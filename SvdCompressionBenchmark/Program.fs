namespace SvdCompressionBenchmark

module Program =

    open MathNet.Numerics
    open MathNet.Numerics.LinearAlgebra
    open RedSvdDriver

    let sqr x = x * x

    let calculate f (input : Matrix<double>) =
        (* SVD *)
        let stopWatch = System.Diagnostics.Stopwatch.StartNew()
        let singularValues, u, vt = f input // run calculation
        stopWatch.Stop()
        //let reconstructed : Matrix<double> = u * singularValues * vt
        //let MSE = Distance.MSE(input.ToColumnWiseArray(), reconstructed.ToColumnWiseArray())
        //let range = (input.Enumerate() |> Seq.max) - (input.Enumerate() |> Seq.min)
        //let NRMSD = (sqrt MSE) / range
        //let PSNR = 20.0 * log10 (1.0 / NRMSD)
        //printfn "  MSE: %A" MSE
        //printfn "  NRMSD: %A" NRMSD
        //printfn "  PSNR: %A" PSNR
        printfn "  Execution time: %A" stopWatch.Elapsed

    let trim rank (s : Matrix<double>, u : Matrix<double>, vt : Matrix<double>) =
        match rank with
        | 0 -> (CreateMatrix.Diagonal(1, 1, 0.0), CreateMatrix.Dense(u.RowCount, 1, 0.0), CreateMatrix.Dense(1, vt.ColumnCount, 0.0))
        | _ -> let trimmedS = s.Diagonal() |> Seq.take rank |> Seq.toArray |> CreateMatrix.Diagonal
               let trimmedU = CreateMatrix.DenseOfColumnVectors (u.EnumerateColumns() |> Seq.take rank)
               let trimmedVT = CreateMatrix.DenseOfRowVectors (vt.EnumerateRows() |> Seq.take rank)
               (trimmedS, trimmedU, trimmedVT)

    let wowFeature () =
        // random matrix with standard distribution:
        let A = (DenseMatrix.randomStandard<double> 100 100)
        let s, u, vt = RedSvdDriver.computeSvdExact A
        let alpha = A.ToColumnWiseArray() |> Array.map sqr |> Array.sum
        let beta = s.Diagonal().Enumerate() |> Seq.map sqr |> Seq.sum
        printfn "alpha: %f" alpha
        printfn "beta: %f" beta // alpha = beta !!! wow

    let benchmark count rows columns =
        for i in 1 .. count do
            printfn "=== Iteration %i ===" i
            printfn "EXACT"
            let randomInput = DenseMatrix.randomStandard<double> rows columns
            calculate computeSvdExact randomInput
            printfn "RANDOMIZED"
            let rank = min rows columns
            calculate (computeSvdRandomized rank) randomInput


    [<EntryPoint>]
    let main argv =

       benchmark 1 100 10000

       System.Console.ReadLine() |> ignore

       0 // return an integer exit code
