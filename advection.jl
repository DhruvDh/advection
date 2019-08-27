function advection(N, dt)
    tmax = 2.0
    nbstep = (tmax / dt)
    xmin = 0.0
    xmax = 1.0
    dx = (xmax - xmin) / 2.0
        
    v = 1.0
    xc = 0.25
    
    α = v * dt/(2.0 * dx)

    x = map(i -> xmin + (i-1)*dx, 1:(N+3))
    u0 = map(i -> -200 * (x[i] - xc)^2, 1:(N+3))

    u = copy(u0)
    unew = copy(u0)

    for timestamp = 1:nbstep
        current_time = timestamp * dt

        for j = 2:(N+2)
            unew[j] = u[j] - α*(u[j+1] - u[j-1]) + 0.5*(u[j+1] - 2*u[j] + u[j-1])
        end
        
        u = copy(unew)
        u[1] = u[N+2]
        u[N+3] = u[2]
    end
end

using BenchmarkTools

case_1_time = @belapsed advection(103, 0.0009)
case_2_time = @belapsed advection(1003, 0.00009)

println("Case 1: N = 103, dt = 0.0009 took $case_1_time seconds")
println("Case 2: N = 1003, dt = 0.00009 took $case_2_time seconds")