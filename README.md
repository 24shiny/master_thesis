# BRIEF INTRODUCTION
&nbsp;&nbsp; These are some selected MATLAB codes for my master thesis ***Extrapolation of Quantum TIme Series***. Its full-text is available on <strong>u:thesis</strong>, a ditigal archive run by <strong>the University of Vienna</strong>:<br>
https://utheses.univie.ac.at/search/?authors=shin&title_subtitle=shin&strict_search=false&sortResultsField=score&sortResultsOrder=desc&resultsPerPage=25&currentPage=1
<br>

# ABSTRACT
<p align="justify">
&nbsp;&nbsp; It is one of the primary goals of physics to predict what a physical quantity of interest is going to be like in the future based on currently available information. In this sense, extrapolation is along the lines of physics. In the field of numerical analysis, extrapolation refers to a method to extend a given set of data points and estimate a value beyond the range. In this study, narrowing down our focus to Hilbert-Schmidt (HS) observables, we try to extrapolate their time averages. To this end, we come up with two possible scenarios and two ‘HS extrapolation functions’, based on superoscillations and McLaurin expansion, respectively. While trying to minimize the error associated with the extrapolation, we find that both scenarios boil down to a convex optimization over extrapolation functions with minimum l1 norm. For the sake of simplicity, a functional ansatz in the form of a series is adopted and the optimization problem is solved with respect to its coefficient vector by utilizing the software packages named MOSEK and CVX. Remarkably, optimal extrapolation functions feature ‘sparse coefficients’, namely, only few non-zero components. In this regard, we first study how their indices change with the estimation error δ, which allows us to group them. At the same time, from the fact that the number of sparse coefficients increases as δ decreases, we deduce that the coefficient vector ends up with mostly non-zero components as δ approximates zero. On the other hand, linear fit functions to values of each group reveal that their values can be predicted with high accuracy. Lastly, we extend the ansatz to a future time point  τ, by repeating the same steps with respect to τ. During this process, a new index group is discovered, which strengthens our confidence in the deduction. Afterwards, we turn our attention to the error model, which indicates how reliable each point in the time series is. In order to figure out the correspondence between extrapolation function and error models, we reformulate the original optimization problem based on its duality and consider three extrapolation functions: Lagrange polynomials and the two HS extrapolation functions mentioned above. The recast problem takes as its objective function the duality gap, which serves as an indicator of the optimality of a given function with the error model. We thereby compare the optimality of the candidate functions based on their optimal values. It turns out that none of the considered functions is optimal for any error model. However, they are all close to optimal for the zero-error model.
<br><br>
<strong>Keywords</strong> : Numerical analysis, Extrapolation, Matrix state approximation, Hilbert-Schmidt observables, Convex optimization, CVX, MOSEK. 
</p>

# CORRECTION 
&nbsp;&nbsp; Some typos were belately found:
+ li
+ li
+ li

To be updated! 😎
