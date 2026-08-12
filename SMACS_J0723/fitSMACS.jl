using LensFactory

# Parse the YAML file
parameters = LensModel.read_input("SMACS0723_input.yaml")

time_start = time()
LensModel.fit_model(parameters)
time_end = time()

println("Total time taken: $(time_end - time_start) seconds")
