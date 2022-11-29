

## Coding Conventions

This software is organised into various modules.
Each module has various type definitions and methods declared.
Commenting is 
done inside each module file to describe the purpose of the module and its usage.
Most of the comments in the code are inline with the text in this documentation. 
Within each module, variable naming is done 
to reduce the effort needed to understand the source code.
For example, 
`distance = velocity * time` is prefered over using 
`a = b * c` in most parts of the software.
The code inside each method is properly intended using spaces to facilitate 
redability.
We followed this documentation [![guide](http://www.stochasticlifestyle.com/finalizing-julia-package-documentation-testing-coverage-publishing/)](http://www.stochasticlifestyle.com/finalizing-julia-package-documentation-testing-coverage-publishing/)


* The methods ending with `!` ideally should not allocate and memory. They are supposed to be fast and iteratively called inside loops.


