namespace SvdCompressionBenchmark

module RedSvdDriver =

    module Native =
    
        open System.Runtime.InteropServices
        open Microsoft.FSharp.NativeInterop
    
        [<DllImport("redsvd.dll")>]
        extern void ComputeSvdExact(double[] inputMatrix_RowMajor, int numberOfRows, int numberOfColumns, [<Out>] double[] singularValues, [<Out>] double[] U_VT_ColumnMajor)

        [<DllImport("redsvd.dll")>]
        extern void ComputeSvdRandomized(double[] inputMatrix_RowMajor, int numberOfRows, int numberOfColumns, int rank, [<Out>] double[] singularValues, [<Out>] double[] U_VT_ColumnMajor)
    
    open MathNet.Numerics.LinearAlgebra

    let computeSvdExact (A : Matrix<double>) = 
        
        let rank = min A.RowCount A.ColumnCount
        
        let singularValues = Array.zeroCreate rank
        let U_VT_ColumnMajor = Array.zeroCreate (A.RowCount * rank + rank * A.ColumnCount)
        
        Native.ComputeSvdExact(A.ToRowWiseArray(), A.RowCount, A.ColumnCount, singularValues, U_VT_ColumnMajor)
        
        (CreateMatrix.Diagonal singularValues, CreateMatrix.Dense (A.RowCount, rank, U_VT_ColumnMajor.[.. (A.RowCount * rank - 1)]), CreateMatrix.Dense (rank, A.ColumnCount, U_VT_ColumnMajor.[(A.RowCount * rank) ..]))


    let computeSvdRandomized rank (A : Matrix<double>) = 
        
        let singularValues = Array.zeroCreate rank
        let U_VT_ColumnMajor = Array.zeroCreate (A.RowCount * rank + rank * A.ColumnCount)
        
        Native.ComputeSvdRandomized(A.ToRowWiseArray(), A.RowCount, A.ColumnCount, rank, singularValues, U_VT_ColumnMajor)

        (CreateMatrix.Diagonal singularValues, CreateMatrix.Dense (A.RowCount, rank, U_VT_ColumnMajor.[.. (A.RowCount * rank - 1)]), CreateMatrix.Dense (rank, A.ColumnCount, U_VT_ColumnMajor.[(A.RowCount * rank) ..]))

