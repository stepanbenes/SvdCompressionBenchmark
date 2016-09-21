namespace SvdCompressionBenchmark

module Program =

    open MathNet.Numerics.LinearAlgebra
    open RedSvdDriver

    let calculate f input =
        (* SVD *)
        let singularValues, u, vt = f input
        
        printfn "Singular Value Decomposition of:"
        printfn "%A" input
        
        printfn "singular values:"
        printfn "%A" singularValues
        printfn "U:"
        printfn "%A" u
        printfn "VT:"
        printfn "%A" vt

        let assembled = u * singularValues * vt

        printfn "assembled matrix:"
        printfn "%A" assembled

    let benchmark count =
        
        for i in 1 .. count do
            let randomInput = DenseMatrix.randomStandard<double> 10 10
            calculate computeSvdExact randomInput
            calculate (computeSvdRandomized 9) randomInput

    [<EntryPoint>]
    let main argv = 
        
        benchmark (System.Int32.Parse argv.[0])

        // random matrix with standard distribution:
        //let m6 = DenseMatrix.randomStandard<double> 3 4
        
//        let input = matrix [[ 1.0;  -2.0;   3.0];
//                            [ 5.0;   8.0;  -1.0];
//                            [ 2.0;   1.0;   1.0];
//                            [-1.0;   4.0;  -3.0]]
        
        

        (* Randomized SVD *)
        //let singularValues, uv, t = RedSvdDriver.computeSvdRandomized input 2

        0 // return an integer exit code

    
