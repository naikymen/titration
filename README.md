
# Numerical method

It solves the acid-base equations for a given Na concentration, in a phosphoric acid buffer.
 
`R/phosphoric_titration_numeric.R` implements the numerical method, and was run with R version 4.0.4, and `nleqslv` version 3.3.2.

Thanks to Bhas for helping with that at stackoverflow: https://stackoverflow.com/a/66883518/11524079

The result is saved for convenience and reference in the `output/result.df.RDS` file (and others, similarly named).

Everything relevant in a plot :)

![numerical plot](./output/numerical.png)

> Titration curve (top left); values of objective functions (top right, BM: mass balance, BC: charge balance, An: acid-base equilibria); speciation diagram (bottom left) and a nice _Na VS species_ plot which I never saw in class (bottom right).

# Analytical method

`R/phosphoric_titration_analytic.R` implements the analytical method explained by Poutnik at chem.se: https://chemistry.stackexchange.com/a/149285/107836

Very simple and elegant: solving for Na is _much_ simpler than solving for pH.

It took me 10 minutes to implement Poutnik's solution (while learning how to use numerical solvers took me days).

> Moraleja: KISS :P

# Experimental data

Code for comparison to experimental data (from Julia Martín _et al_, DOI 10.20431/2349-0403.0409002) is at `R/comparison.R`.

![comparison plot](./output/comparison.png)

Well, its something!

