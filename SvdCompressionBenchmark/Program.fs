namespace SvdCompressionBenchmark

module Program =

    open MathNet.Numerics
    open MathNet.Numerics.LinearAlgebra
    open RedSvdDriver

    let sqr x = x * x

    let calculate f (input : Matrix<double>) =
        (* SVD *)
        let singularValues, u, vt = f input
        
//        printfn "Singular Value Decomposition of:"
//        printfn "%A" input
//        
//        printfn "singular values:"
//        printfn "%A" singularValues
//        printfn "U:"
//        printfn "%A" u
//        printfn "VT:"
//        printfn "%A" vt

        let reconstructed : Matrix<double> = u * singularValues * vt

//        printfn "assembled matrix:"
//        printfn "%A" reconstructed

        let mse = Distance.MSE(input.ToColumnWiseArray(), reconstructed.ToColumnWiseArray())
        let psnr = 10.0 * log10 ((sqr 255.0) * mse)

        printfn "MSE: %A" mse
        printfn "PSNR: %A" psnr

    let benchmark count =
        
        for i in 1 .. count do
            let randomInput = DenseMatrix.randomStandard<double> 100 100
            calculate computeSvdExact randomInput
            calculate (computeSvdRandomized 100) randomInput

    let trim rank (s : Matrix<double>, u : Matrix<double>, vt : Matrix<double>) =
        match rank with
        | 0 -> (CreateMatrix.Diagonal(1, 1, 0.0), CreateMatrix.Dense(u.RowCount, 1, 0.0), CreateMatrix.Dense(1, vt.ColumnCount, 0.0))
        | _ -> let trimmedS = s.Diagonal() |> Seq.take rank |> Seq.toArray |> CreateMatrix.Diagonal
               let trimmedU = CreateMatrix.DenseOfColumnVectors (u.EnumerateColumns() |> Seq.take rank)
               let trimmedVT = CreateMatrix.DenseOfRowVectors (vt.EnumerateRows() |> Seq.take rank)
               (trimmedS, trimmedU, trimmedVT)

    [<EntryPoint>]
    let main argv = 
        
        //benchmark 100

       let size = 100
       
       // random matrix with standard distribution:
       let A = (DenseMatrix.randomStandard<double> size size)
       
       //let A = matrix [[ 1.0;  -2.0;   3.0];
       //                    [ 5.0;   8.0;  -1.0];
       //                    [ 2.0;   1.0;   1.0];
       //                    [-1.0;   4.0;  -3.0]]
       
       let s, u, vt = RedSvdDriver.computeSvdExact A
       
       let alpha = A.ToColumnWiseArray() |> Array.map sqr |> Array.sum
       let beta = s.Diagonal().Enumerate() |> Seq.map sqr |> Seq.sum
       printfn "alpha: %f" alpha
       printfn "beta: %f" beta // alpha = beta !!! wow

//       for rank in 0 .. size do
//
//            printfn ">> rank = %A" rank
//            
//            let s2, u2, vt2 = trim rank (s, u, vt)
//
//            let B = u2 * s2 * vt2
//            let D = A - B
//
//            //printfn "D = %A" A
//
//            let mn = float(A.RowCount) * float(A.ColumnCount)
//
//            let mseD = (D.ToColumnWiseArray() |> Array.map sqr |> Array.sum) / mn
//            let mseS = (s.Diagonal().Enumerate() |> Seq.skip rank |> Seq.map sqr |> Seq.sum) / mn
//
//            printfn "MSE(D): %.14f" mseD
//            printfn "MSE(S): %.14f" mseS
//            
//            let RMSD = sqrt mseD
//
//            printfn "RMSD: %.14f" RMSD
//
//            let range = (B.Enumerate() |> Seq.max) - (B.Enumerate() |> Seq.min)
//            let NRMSD = RMSD / range
//
//            let PSNR = 20.0 * (log10 0.1)
//
//            printfn "NRMSD: %.14f" NRMSD
//            printfn "PSNR: %.14f" PSNR
//
       System.Console.ReadLine() |> ignore

       0 // return an integer exit code
