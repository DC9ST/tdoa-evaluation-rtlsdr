function [ output_real ] = remove_mean( input_real )
    output_real = input_real - mean(input_real);
end

