# EigenSparseProduct
Made it multi-threading: C = A * B, where A is an `Eigen::ColMajor` sparsematrix and B is an `Eigen::RowMajor` sparsematrix. <br />
实现了`Eigen::ColMajor`稀疏矩阵乘`Eigen::RowMajor`稀疏矩阵的多线程乘法。

### Benckmark
Here we made 2 test cases: 
* A is `471538`\*`941014`, B is `941014`\*`471538`, with `0.00127%` nonzero values.
* A is `235771`\*`235771`, B is `235771`\*`235771`, with `0.00296%` nonzero values.

The machine has an Intel Core-i7 6700K CPU @ 4.38GHz, with 16GiB 2400MHz DDR4 memory (single channel) on Windows 10.0.17134.407. We made the product of each case continuously for 4 times, and logged the total time consumed as following.

<table>
		<tr>
			<td>Method</td>
			<td>Case</td>
			<td>1st</td>
			<td>2nd</td>
			<td>3rd</td>
			<td>4th</td>
		</tr>
		<tr>
			<td>Eigen's</td>
			<td>1</td>
			<td>1141ms</td>
			<td>1176ms</td>
			<td>1162ms</td>
			<td>1162ms</td>
		</tr>
		<tr>
			<td><strong>Ours</strong></td>
			<td>1</td>
			<td>345ms</td>
			<td>298ms</td>
			<td>313ms</td>
			<td>298ms</td>
		</tr>
		<tr>
			<td>Eigen's</td>
			<td>2</td>
			<td>471ms</td>
			<td>479ms</td>
			<td>476ms</td>
			<td>475ms</td>
		</tr>
		<tr>
			<td><strong>Ours</strong></td>
			<td>2</td>
			<td>140ms</td>
			<td>115ms</td>
			<td>121ms</td>
			<td>120ms</td>
		</tr>
</table>

This method turns out to perform a little worse for the 1st product due to memory allocation, but it becomes faster for the following. <b>The overall speed-up ratio is 3.30 ~ 4.16 on this 4 cored machine, according to the table.</b>
