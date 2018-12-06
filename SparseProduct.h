#pragma once
#include <Eigen\Sparse>
#include <vector>
#include <mutex>

//实现并行快速稀疏矩阵乘法C=A*B
//矩阵A需要为列优先，矩阵B需要为行优先
//https://www.licc.vip/article.aspx?id=63
template<typename Scalar, int IsRowMajor = Eigen::SparseMatrix<Scalar>::IsRowMajor>
class SparseProduct
{
private:
	std::vector<std::shared_ptr<std::mutex>> BucketMutex;
	std::vector<std::vector<Eigen::Triplet<Scalar>>> BucketData;
	Eigen::VectorXi Sizes;
	int Rows, Cols;
	int Bucket;
	int Stride;

public:
	SparseProduct(int rowsA, int colsB, int bucketSize = 8192) :
		Rows(rowsA), Cols(colsB), Bucket(bucketSize), Stride(IsRowMajor ? rowsA / Bucket + 1 : colsB / Bucket + 1)
	{
		BucketMutex.resize(bucketSize);
		BucketData.resize(bucketSize);
		if (MySparseMatrix::IsRowMajor) {
			Sizes.resize(rowsA);
		}
		else {
			Sizes.resize(colsB);
		}

		for (int i = 0; i < bucketSize; i++) {
			BucketMutex[i] = std::shared_ptr<std::mutex>(new std::mutex());
		}
	}

	virtual ~SparseProduct()
	{
	}

public:
	//Reserve estimated spaces for BucketData 
	void Reserve(int mid, int nonZerosA, int nonZerosB) {
		for (int i = 0; i < Bucket; i++) {
			int estimatedCountPerLine = ((long long)nonZerosA) * ((long long)nonZerosB) / mid;
			int capacity = estimatedCountPerLine * 2 / Bucket;

			if (BucketData[i].capacity() < capacity)
				BucketData[i].reserve(capacity);
		}
	}

	//Calculate C = A * B
	Eigen::SparseMatrix<Scalar, IsRowMajor> & Product(Eigen::SparseMatrix<Scalar, IsRowMajor> & C, const Eigen::SparseMatrix<Scalar, Eigen::ColMajor> & A, const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> & B) {
		const int rows = this->Rows;
		const int cols = this->Cols;
		const int stride = this->Stride;
		const int len = A.cols();
		if (A.rows() != rows) return C;
		if (B.cols() != cols) return C;
		if (B.rows() != len) return C;
		C.resize(rows, cols);

		Reserve(len, A.nonZeros(), B.nonZeros());
		for (int i = 0; i < BucketData.size(); i++) BucketData[i].clear();

#pragma omp parallel
		{
			//Insert non-zero components to BucketData
#pragma omp for schedule(dynamic, 100)
			for (int i = 0; i < len; i++) {
				if (IsRowMajor) {
					for (MyColSparseMatrix::InnerIterator it1(A, i); it1; ++it1) {
						const int row = it1.row();
						const int bucket = row / stride;
						BucketMutex[bucket]->lock();
						for (MyRowSparseMatrix::InnerIterator it2(B, i); it2; ++it2) {
							const int col = it2.col();
							auto t = Eigen::Triplet<MyMesh::Scalar>(row, col, it1.value() * it2.value());
							BucketData[bucket].push_back(t);
						}
						BucketMutex[bucket]->unlock();
					}
				}
				else {
					for (MyRowSparseMatrix::InnerIterator it2(B, i); it2; ++it2) {
						const int col = it2.col();
						const int bucket = col / stride;
						BucketMutex[bucket]->lock();
						for (MyColSparseMatrix::InnerIterator it1(A, i); it1; ++it1) {
							const int row = it1.row();
							auto t = Eigen::Triplet<MyMesh::Scalar>(row, col, it1.value() * it2.value());
							BucketData[bucket].push_back(t);
						}
						BucketMutex[bucket]->unlock();
					}
				}
			}

			//Reserve C according to BucketData
			std::vector<std::vector<int>> sets(stride);			//using vector instead of unordered_map for detecting duplication turns out be faster, for elements with no more than 50
#pragma omp single
			{
				Sizes.setZero();
			}
#pragma omp for schedule(dynamic, 10)
			for (int b = 0; b < Bucket; b++) {
				for (int i = 0; i < stride; i++) sets[i].clear();

				for (auto & data : BucketData[b]) {
					const int row = data.row();
					const int col = data.col();

					if (IsRowMajor) {
						bool existed = false;
						auto & set = sets[row % stride];
						for (auto & val : set) {
							if (val == col) {
								existed = true;
								break;
							}
						}
						if (!existed) set.push_back(col);
					}
					else {
						bool existed = false;
						auto & set = sets[col % stride];
						for (auto & val : set) {
							if (val == row) {
								existed = true;
								break;
							}
						}
						if (!existed) set.push_back(row);
					}
				}

				for (int i = 0; i < stride && i + b * stride < Sizes.size(); i++) {
					Sizes[i + b * stride] = sets[i].size();
				}
			}
#pragma omp single
			{
				C.reserve(Sizes);
			}

			//Assign C according to BucketData
#pragma omp for schedule(dynamic, 10)
			for (int b = 0; b < Bucket; b++) {
				for (auto & data : BucketData[b]) {
					C.coeffRef(data.row(), data.col()) += data.value();
				}
			}
		}

		return C;
	}
};
