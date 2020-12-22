This is a basic example of conduction-diffusion equation.

-e*\laplace u + w dotmultiply grad u = f ;

u(x,y) = x*((1-exp((y-1)/e))/(1-exp(-2/e)));

output.cpp is a file to check if the output in the test.cpp is correct!

Steps:
1. easymesh D;
2. ./Command // This is a shell script to help compile the .cpp;
3. ./run D to get the .dx files;
4. dx   Using the OpenDX to have see the result of the calculation;
