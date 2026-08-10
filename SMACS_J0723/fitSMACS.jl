using LensFactory

# Parse the YAML file
parameters = LensModel.read_input("SMACS0723_input.yaml")

LensModel.fit_model(parameters)
