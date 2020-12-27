---
title: "Kernel Density Estimation"
date: "2020-12-26"
output: 
  html_document:
    keep_md: true
---



# Introduction

One of the last assignments I did before graduation was a tutorial on Kernel Density Estimation (KDE) for my classmates taking Mathematical Statistics. This is the tutorial and activities I wrote with my classmate, Kaden Bieger, to walk our class through the basics of KDE. Working through the activities in this tutorial should give you a good working understanding of Kernel Density Estimation!


Nonparametric statistics is a rapidly developing field. It is also very different than the parametric content covered in most Introduction to Mathematical Statistics classes! Broadly speaking, nonparametric methods allow us to relax assumptions about our data. We often assume our data comes from a normal distribution, or at the very least from a distribution with mean $\mu$ and variance $\sigma^2$. Nonparametric methods are not based on such parameters.

*Kernel density estimation* (KDE) is a common technique used to estimate probability density functions (PDFs). In practice, we rarely know much at all about the true distribution of our sampled data. KDE allows us to get an estimated PDF without making unreasonable assumptions of our data. 

## Intuition

KDE is a method you've seen before even if you don't know it! The intuition behind KDE is similar to the intuition behind a simple **histogram**. Histograms can be used to visually approximate the distribution of data. There are a few reasons we want a more sophisticated method, however. First, as the plot below illustrates, histograms are very sensitive to the number and width of bins. Additionally, a histogram provides an excessively local estimate -- there is no "smoothing" between neighboring bins. 



```r
set.seed(455)
simdat <- rnorm(n = 100, mean = 2, sd = 4)
# fn borrowed from Son Phan
hist_bins <- function(data, bins) {
  ggplot(mapping = aes(x = data, y=..density..)) +
  geom_histogram(bins = bins, fill = "steelblue3") + 
  ggtitle(sprintf(fmt = "%i bins", bins)) + theme_minimal()
}

grid.arrange(hist_bins(simdat, 10), 
             hist_bins(simdat, 20), 
             hist_bins(simdat, 40), 
             hist_bins(simdat, 80), nrow=2)
```

![](introduction-to-KDE_files/figure-html/unnamed-chunk-2-1.png)<!-- -->


If you've ever made a **density plot**, you have even more experience with KDE. Behind the scenes, the `ggplot2` function `geom_density` employs KDE! It turns out we can set the specific **kernel** and **bandwidth** when we call `geom_density`. We'll cover kernels and bandwidths below.


# Kernel Density Estimation

The goal of kernel density estimation to to estimate the PDF of our data without knowing its distribution. We use a **kernel** as a weighting function to smooth our data in order to get this estimate. 

A **kernel** is a probability density function with several additional conditions:

1. Kernels are non-negative and real-values
2. Kernels are symmetric about 0

Several familiar PDFs, including the Gaussian and Uniform PDFs, meet these requirements. Once we have a kernel selected, we can implement KDE.

Given data $x = (x_1, x_2, \ldots, x_n)$, a kernel function $K$, and a selected bandwidth $h$, the kernel density estimate at point $s$ is defined

\begin{equation}
  \hat{f}_n(x) = \frac{1}{nh} \sum_{i = 1}^{n} K\left(\frac{x - x_i}{h}\right)
\end{equation}

This is the kernel density estimator at a *single* point. To estimate an entire PDF, we apply the kernel to each point in our sample. The procedure is as follows: 

1. Apply the kernel function each data point in our sample
![](introduction-to-KDE_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

2. Sum all $n$ kernel functions
![](introduction-to-KDE_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

3. Divide by $n$. Because each kernel function integrates to one, the resulting KDE will still integrate to one. 
![](introduction-to-KDE_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

In notation, we can write this procedure as $f(x) = \frac{1}{n} \sum_{i=1}^{n} K(x - x_i)$ where $x-x_i$ centers the kernel function on each $x_i$. Note the similarity to the formal definition of a kernel density estimator above. The only additional parameter is the bandwidth!




#### Acitivity 1

Consider 100 independent and identically distributed samples from an unknown distribution (plotted below). Given a bandwidth $h = 0.5$, write the Gaussian kernel density estimator for this sample at $x_i = 5$. We can use a standard normal with $\sigma^2 = 1$ and $\mu = 0$. 


```r
ggplot(asimdat, aes(x=x)) +
  geom_histogram(bins = 15, fill = "steelblue3") +
  theme_minimal() + ggtitle("iid sample from unknown distribution")
```

![](introduction-to-KDE_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

Try plotting your calculated kernel density estimate below. There are lots of built in ways to calculate KDEs, but we've provided a code outline so you can see how this process works "under the hood." Compare your results to output of `geom_density`!


```r
# delete eval = FALSE above if you want the plots to show up in your knitted doc!

my_est <- function(x, xi = 5, h = 0.5){
  # put your estimate here!
}


# to get the whole kde estimate, sum my_est over each s
kern <- matrix(ncol = 100, nrow = 100)

kde <- function(x, xi = 5, h = 0.5){
  for(i in 1:100){
    kern[i, ] <- my_est(x = ???, s = ???, h = h)
  }
  colSums(kern)
}


# uncomment lines in this ggplot code as you finish the functions above

ggplot(asimdat, aes(x=x)) +
  theme_minimal() +
  geom_density() +   # generic geom_density
  geom_density(bw = 0.5, kernel = "gaussian", color = "blue") + #closer to your estimator
#  stat_function(fun = my_est, color = "red")  + # your estimate at s = 5
#  geom_line(aes(x = x, y = kde(x)), color = "pink") + # and your overall estimate!
  NULL 
```




# Choosing a Kernel Function

Kernel functions are non-negative, symmetric, and decreasing for all $x > 0$ (because they're symmetric, this also means they are increasing for all $x < 0$). Gaussian and Rectangular (uniform) are two kernel choices, though there are many others. The plot below shows kernel density estimates for the same simulated data using four common kernels. In this example, we use a fairly low bandwidth (0.3), to make differences between the kernels clear. 


```r
gaus <- ggplot(mapping = aes(x = simdat, y=..density..)) +
  geom_histogram(bins = 10, fill = "steelblue3") + theme_minimal() +
  geom_density(kernel = "gaussian", bw = 0.3, color = "pink", size = 1) +
  ggtitle("Gaussian Kernel")

e <- ggplot(mapping = aes(x = simdat, y=..density..)) +
  geom_histogram(bins = 10, fill = "steelblue3") + theme_minimal() +
  geom_density(kernel = "epanechnikov", bw = 0.3, color = "pink", size = 1) +
  ggtitle("Epanechnikov Kernel")

tri <- ggplot(mapping = aes(x = simdat, y=..density..)) +
  geom_histogram(bins = 10, fill = "steelblue3") + theme_minimal() +
  geom_density(kernel = "triangular", bw = 0.3, color = "pink", size = 1) +
  ggtitle("Triangular Kernel")

rect <- ggplot(mapping = aes(x = simdat, y=..density..)) +
  geom_histogram(bins = 10, fill = "steelblue3") + theme_minimal() +
  geom_density(kernel = "rectangular", bw = 0.3, color = "pink", size = 1) +
  ggtitle("Rectangular Kernel")

grid.arrange(gaus, e, tri, rect, nrow=2)
```

![](introduction-to-KDE_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

Choice of kernel doesn't affect end result that much, although there are benefits to some kernels in particular. The Epanechnikov kernel, for example, is MSE optimal. 

**Bonus Activity:** If you want to try out other kernel options inside `geom_density`, go to "Help" in Rstudio and search for "density." Scroll down to "arguments" to see all possible kernels!


# Choosing Bandwidth

The bandwidth is the width of our smoothing window -- following the histogram example, this is like the width of bins. Bandwidth is generally represented as $h$. Bandwidth is how we control the smoothness of our estimate. The figure below shows a Gaussian KDE with various bandwidths. 


```r
gaus25 <- ggplot(mapping = aes(x = simdat, y=..density..)) +
  geom_histogram(bins = 10, fill = "steelblue3") + theme_minimal() +
  geom_density(kernel = "gaussian", bw = 0.25, color = "pink", size = 1) +
  ggtitle("Gaussian Kernel, bandwith = 0.25")

gaus50 <- ggplot(mapping = aes(x = simdat, y=..density..)) +
  geom_histogram(bins = 10, fill = "steelblue3") + theme_minimal() +
  geom_density(kernel = "gaussian", bw = 0.50, color = "pink", size = 1) +
  ggtitle("Gaussian Kernel, bandwith = 0.50")

gaus75 <- ggplot(mapping = aes(x = simdat, y=..density..)) +
  geom_histogram(bins = 10, fill = "steelblue3") + theme_minimal() +
  geom_density(kernel = "gaussian", bw = 0.75, color = "pink", size = 1) +
  ggtitle("Gaussian Kernel, bandwith = 0.75")

gaus1 <- ggplot(mapping = aes(x = simdat, y=..density..)) +
  geom_histogram(bins = 10, fill = "steelblue3") + theme_minimal() +
  geom_density(kernel = "gaussian", bw = 1, color = "pink", size = 1) +
  ggtitle("Gaussian Kernel, bandwith = 1")

grid.arrange(gaus25, gaus50, gaus75, gaus1, nrow=2)
```

![](introduction-to-KDE_files/figure-html/unnamed-chunk-10-1.png)<!-- -->


## Bandwidth tradeoffs

A small $h$ places more weight on the data, while a large $h$ performs more smoothing. We can use smaller bandwidths when $n$ is larger, particularly when the data has a small range (ie it is tightly packed). Larger $h$ is better when we have less data or more spread out observations. 

There's an analogy here with Bayesian priors -- choosing $h$ lets us choose how much weight we want to give the data. 


#### Activity 2

Using the kernel you defined above, plot a KDE estimate using 3 different values for $h$. You can do so by updating the bandwidth you defined in `my_est`.


```r
# uncomment this ggplot code to get started

#ggplot(asimdat, aes(x = x, y = kde(x))) + # play with the kde function! h = ?
#  geom_line() + theme_minimal()
```



# Bias, Variance, MSE

## Bias
Naturally, we would like to consider the bias and mean squared error of estimators produced by kernel density estimation. In this section, we outline proofs deriving the expected value, bias, and mean squared error for kernel density estimates. The proofs here are somewhat simplified but largely rely on **u-substitution** and **Taylor expansions**. 

Given observations $x_1, x_2, \ldots, x_n \overset{iid}{\sim} f$, we define the expected value of our estimator as follows: 

$$
\mathbf{E}[\hat{f}_n(x)] = \mathbf{E}\left[\frac{1}{nh} \displaystyle \sum_{i = 1}^{n} K \left( \frac{x - x_i}{h}\right) \right] = \frac{1}{h}\mathbf{E}\left[K \left( \frac{x - X}{h}\right) \right] = \frac{1}{h} \int K\left( \frac{x - X}{h}\right) f(x) dx
$$

This first line follows simply from the definitions of kernel density estimators and expected value. Here, we can let $w = X$ so we have

$$
\mathbf{E}[\hat{f}_n(x)] = \frac{1}{h} \int K\left( \frac{x - w}{h}\right) f(w) dw
$$

Next, we let $u = \frac{x - w}{h}$ and substitute. Note we can rearrange $u = \frac{x - w}{h}$ to get $w = x - hu$.

$$
\mathbf{E}[\hat{f}_n(x)] = \frac{1}{h} \int K\left( u\right) f(x - hu) du
$$

Then, we apply the 2nd order Taylor expansion for $f(x - hu)$ about $h = 0$. Omitting several steps of algebra, our Taylor expansion can be written

$$
\begin{aligned}
f(x - hu) &=  f(x) + \frac{f'(x)}{1!}(u)(h-0) + \frac{f''(x)}{2!}(u^2)(h-0)^2 + o(h^2) \\
    &= f(x) - huf'(x) + \frac{h^2u^2}{2}f''(x) + o(h^2)
\end{aligned}
$$

where $o(h^2)$ is some function which approaches zero *more quickly* than $h^2$ as $h^2$ approaches 0. Technically, this $o()$ notation indications that $o(h^2)$ becomes negligible compared to $h^2$ as $h \rightarrow 0$. **For our purposes, we assume $o(h^2) \rightarrow 0$ as $h^2 \rightarrow 0$**.

We plug the Taylor expansion into our expected value above and simplify via algebra. 

$$
\begin{aligned}
\mathbf{E}[\hat{f}_n(x)] & = \int K(u) \left[f(x) - huf'(x) + \frac{h^2u^2}{2}f''(x) + o(h^2)\right] du \\ \\
& = f(x)\int K(u) du + hf'(x)\int uK(u) du +  \frac{h^2}{2}f''(x)\int u^2 K(u)du  + o(h^2) \\
& = f(x) + \frac{h^2}{2}f''(x)\int u^2 K(u)du + o(h^2)
\end{aligned}
$$

We can plug this into the definition of bias such that 

$$
\begin{aligned}
\textbf{Bias}(\hat{f}_n(x)) &= E[\hat{f}_n(x)] - f(x) \\
&= \frac{h^2}{2}f''(x)\int u^2 K(u)du + o(h^2) \\
    &= \frac{t \cdot h^2}{2}f''(x) + o(h^2)
\end{aligned}
$$

where $t = \int u^2 K(u)du$. 

#### Activity 3
What happens to bias as bandwidth $h$ increases/decreases? 

## Variance
The proof for variance is very similar to the proof for bias. This proof is possible because $K$ is symmetric about 0. Feel free to work through this proof more on your own!

$$
\begin{aligned}
    \mathbf{Var}(\hat{f}_n(x))
    &=
    \mathbf{Var} \left(\frac{1}{nh} \displaystyle \sum_{i = 1}^{n} K \left(\frac{x - x_i}{h}\right)\right) \\
    &= \frac{1}{nh^2} \left(\mathbf{E}\left[ K^2 \left( \frac{x - x_i}{h}\right) \right] - \mathbf{E}\left[ K \left( \frac{x - x_i}{h}\right)\right]^2\right)\\
    &\leq
    \frac{1}{nh^2} \mathbf{E}\left[ K^2 \left( \frac{x - X}{h}\right) \right] \\
    &=
    \frac{1}{nh^2} \int K^2 \left( \frac{x - w}{h}\right)f(w) dw
  \end{aligned}
$$

As above, we substitute $u = \frac{x-t}{h}$ and plug in a 1st order Taylor expansion of $f(x - hu)$. 

$$
\begin{aligned}
    \mathbf{Var}(\hat{f}_n(x))
    &\leq
    \frac{1}{nh^2} \int K^2(u)f(x - hu)hdu \\
    &=
    \frac{1}{nh} \int K^2(u)f(x - hu)du \\
    &=
    \frac{1}{nh} \int K^2(u)[f(x) + huf'(x) + o(h)]du  \\
    &=
    \frac{1}{nh} \bigg(f(x)\int K^2(u) du + hf'(x)\int uK^2(u) du + o(h)\bigg) \\
    \mathbf{Var}(\hat{f}_n(x)) &\leq \frac{f(x)}{nh}\int K^2(u) du + o\bigg(\frac{1}{nh}\bigg) \\
    &=
    \frac{z}{nh}f(x) + o\bigg(\frac{1}{nh}\bigg)
  \end{aligned}
$$

where $z = \int K^2(u) du$.

#### Activity 4
What happens to variance as $h$ changes? As $n$ changes?

#### Activity 5
Given the Bias and Variance above, find the Mean Squared Error of our estimator. Recall that the formula for MSE is $Var() + Bias^2()$. 


## Asymptotic MSE

Given the MSE we derived above, we can see that the Asymptotic MSE (AMSE) is $\frac{t^2h^4}{4}\left[f''(x)\right]^2 +  \frac{z}{nh}f(x)$. Often, we use the AMSE to optimize $h$ at a given point. 

#### Activity 6

Try optimizing the AMSE in terms of $h$. **Bonus:** Try finding the optimal $h$ for the sample in Activity 1. 


## An Open Question: AMISE

The optimization in Activity 6 optimizes $h$ at a given $x$ value. A common method for optimizing $h$ across an entire distribution is using the asymptotic mean integrated square error (AMISE). As the name suggests, we can define the AMISE as follows: 

$$
\textbf{AMISE}(\hat{f}_n(x)) = \int \left(\frac{t^2h^4}{4}\left[f''(x)\right]^2 +  \frac{z}{nh}f(x) \right) dx
$$

It's possible to optimize AMISE in terms of $h$ to get the optimal bandwidth across an entire sample. We won't ask you to do this -- it's definitely a challenge problem! In the end, you'd find that the optimal $h$ is dependent upon $\int \left[f''(x)\right]^2 dx$, or the curvature of the underlying PDF. Many recent advancements in KDE studies have been in the realm of estimating AMISE or the underlying curvature in order to optimize $h$. 
