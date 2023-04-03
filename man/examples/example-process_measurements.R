data("alligator")

# Process dataset; vertebra index in "Vertebra" column
alligator_data <- process_measurements(alligator,
                                       pos = "Vertebra")
