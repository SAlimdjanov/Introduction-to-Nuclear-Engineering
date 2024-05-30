# Introduction to Nuclear Engineering

This repository contains solutions to examples and problems in each chapter of <i>Introduction to Nuclear Engineering, 3rd Ed.</i> by John R. Lamarsh and Anthony J. Barratta.

## Contents

Solutions are generated and written to text files based on whether they are example problems or chapter problems.

In the `helpers` directory, conversion factors, relevant data and constants, as well as all relevant formulae are in written in their respective files.

Example and Chapter Problem scripts are in the root directory for the following chapters. The first chapter is omitted, as it is an introductory chapter with no problems or examples present.

2. Atomic and Nuclear Physics

**Notes**: The return type of answers in the solution scripts are as follows:

-   For questions with one answer: `Tuple[Union[str, float], str, int]` A tuple that contains either a string answer (i.e. answers with words. In this case, the next two fields are omitted) or a floating point value accompanied with its units and number of significant figures
-   For multi-part answers (i.e. a), b), c)), a list of these tuples are returned
-   The two answer types above are abstracted in custom types named `SingleAnswer` and `MultiPartAnswer`, respectively.
-   Answers stating `Derivation Question` do not contain numeric solutions. These are symbolic derivation questions that involve pen and paper work

### References

Lamarsh, J. R., & Baratta, A. J. (2001). Introduction to Nuclear Engineering, 3rd Ed. Upper Saddle River, NJ: Prentice Hall.
