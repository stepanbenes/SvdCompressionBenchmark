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

    [<EntryPoint>]
    let main argv = 
        
        benchmark 100

        // random matrix with standard distribution:
        //let m6 = DenseMatrix.randomStandard<double> 3 4
        
        let input = matrix [[ 1.0;  -2.0;   3.0];
                            [ 5.0;   8.0;  -1.0];
                            [ 2.0;   1.0;   1.0];
                            [-1.0;   4.0;  -3.0]]
        
        

        let singularValues, uv, t = RedSvdDriver.computeSvdExact input

        let sqr x = x * x

        let mn = float(input.RowCount) * float(input.ColumnCount)

        let mse1 = (input.ToColumnWiseArray() |> Array.map sqr |> Array.sum) / mn
        let mse2 = (singularValues.Diagonal().Enumerate() |> Seq.map sqr |> Seq.sum) / mn

        printfn "MSE1: %A" mse1
        printfn "MSE2: %A" mse2
        printfn "Equals: %A" (mse1 = mse2)

        System.Console.ReadLine() |> ignore

        0 // return an integer exit code

    
