namespace SvdCompressionBenchmark

module Program =

    open MathNet.Numerics
    open MathNet.Numerics.LinearAlgebra
    open RedSvdDriver

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
        let psnr = 10.0 * log10 (255.0**2.0 * mse)

        printfn "MSE: %A" mse
        printfn "PSNR: %A" psnr

    let benchmark count =
        
        for i in 1 .. count do
            let randomInput = DenseMatrix.randomStandard<double> 100 100
            calculate computeSvdExact randomInput
            calculate (computeSvdRandomized 100) randomInput

    let sqr x = x * x

    let trim rank (s : Matrix<double>, u : Matrix<double>, vt : Matrix<double>) =
        let trimmedS = s.Diagonal() |> Seq.take rank |> Seq.toArray |> CreateMatrix.Diagonal
        let trimmedU = CreateMatrix.DenseOfColumnVectors (u.EnumerateColumns() |> Seq.take rank)
        let trimmedVT = CreateMatrix.DenseOfRowVectors (vt.EnumerateRows() |> Seq.take rank)
        (trimmedS, trimmedU, trimmedVT)

    [<EntryPoint>]
    let main argv = 
        
        //benchmark 100

        for i in 1 .. 10 do

            let size = 100
            // random matrix with standard distribution:
            let A = DenseMatrix.randomStandard<double> size size
            
            //let A = matrix [[ 1.0;  -2.0;   3.0];
            //                    [ 5.0;   8.0;  -1.0];
            //                    [ 2.0;   1.0;   1.0];
            //                    [-1.0;   4.0;  -3.0]]
            
            let rank = 70
            let s, u, vt = RedSvdDriver.computeSvdExact A
            let s2, u2, vt2 = trim rank (s, u, vt)

            //let B = u * s * vt
            let B2 = u2 * s2 * vt2

            //let D = A - B
            let D2 = A - B2

            printfn "D2 = %A" D2

            let mn = float(A.RowCount) * float(A.ColumnCount)

            let mseD2 = (D2.ToColumnWiseArray() |> Array.map sqr |> Array.sum) / mn
            let mse2 = (s.Diagonal().Enumerate() |> Seq.skip rank |> Seq.map sqr |> Seq.sum) / mn

            printfn "MSE(D2): %.14f" mseD2
            printfn "MSE2:    %.14f" mse2

            System.Console.ReadLine() |> ignore

        0 // return an integer exit code

    
