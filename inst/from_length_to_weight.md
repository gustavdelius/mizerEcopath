---
title: "Conversion between weight and length"
output:
  pdf_document: default
  html_document: default
---
In mizer we have the following PDE for a density $N(w)$ in a variable w:

$$\partial N / \partial t = - \partial J_w/\partial w - m w^(n-1) N$$

where

$$J_w = (A w^n - B w) N - D w^{n+1} \partial N/\partial w .$$

We want to convert this into a PDE for a density $u(l,t)$ in a variable $l$ related to $w$ by $w= a l ^ b$ where $n = 1 - 1/b$. Write the PDE in the form

$$\partial u / \partial t = - \partial J_l/\partial l - m / l u$$

where

$$J _l= k (L_\infty^0 - l) u - \alpha l \partial u/\partial l .$$

We would like to express the parameters $k$, $L_\infty^0$ and $\alpha$ in terms of the parameters $A$, $B$ and $D$.

### **Step 1: Relate the Densities using Conservation**

The number of particles must be conserved under the change of variables. This means that the number of particles in an infinitesimal interval $dw$ must be equal to the number of particles in the corresponding interval $dl$.
$$N(w,t) dw = u(l,t) dl$$
This gives the relationship between the two densities:
$$u(l,t) = N(w,t) \frac{dw}{dl}$$
Given the transformation $w = a l^b$, we can calculate the derivative $\frac{dw}{dl}$:
$$\frac{dw}{dl} = a b l^{b-1}$$
Substituting this into the relation between the densities, we get:
$$u(l,t) = N(w(l),t) \cdot (a b l^{b-1})$$
From this, we can express $N(w,t)$ in terms of $u(l,t)$:
$$N(w,t) = \frac{u(l,t)}{a b l^{b-1}}$$

### **Step 2: Transform the PDE**

The original PDE is:
$$\frac{\partial N}{\partial t} = - \frac{\partial J_w}{\partial w} - m w^{n-1} N$$
We can write this in conservation form:
$$\frac{\partial N}{\partial t} + \frac{\partial J_w}{\partial w} = - m w^{n-1} N$$
To change variables from $(w,t)$ to $(l,t)$, we multiply the entire equation by $\frac{dw}{dl}$:
$$\frac{\partial N}{\partial t}\frac{dw}{dl} + \frac{\partial J_w}{\partial w}\frac{dw}{dl} = - m w^{n-1} N \frac{dw}{dl}$$
Let's transform each term:

* **Time derivative term:**
    $$\frac{\partial N}{\partial t}\frac{dw}{dl} = \frac{\partial}{\partial t}\left(N \frac{dw}{dl}\right) - N \frac{\partial}{\partial t}\left(\frac{dw}{dl}\right)$$
    Since $\frac{dw}{dl} = a b l^{b-1}$ does not depend on time $t$, the second term is zero. Thus:
    $$\frac{\partial N}{\partial t}\frac{dw}{dl} = \frac{\partial}{\partial t}\left(N \frac{dw}{dl}\right) = \frac{\partial u}{\partial t}$$

* **Flux derivative term:**
    Using the chain rule:
    $$\frac{\partial J_w}{\partial w}\frac{dw}{dl} = \frac{\partial J_w}{\partial l}$$

* **Source term:**
    $$ - m w^{n-1} N \frac{dw}{dl} = - m w^{n-1} u$$
    We use the given relations $w=al^b$ and $n=1-1/b$, which implies $n-1 = -1/b$.
    $$w^{n-1} = (a l^b)^{n-1} = a^{n-1} l^{b(n-1)} = a^{-1/b} l^{-1}$$
    So the source term becomes:
    $$- m a^{-1/b} l^{-1} u = - \frac{m a^{-1/b}}{l} u$$

Substituting these transformed terms back into the PDE, we get:
$$\frac{\partial u}{\partial t} + \frac{\partial J_w}{\partial l} = - \frac{m a^{-1/b}}{l} u$$
Rearranging this gives:
$$\frac{\partial u}{\partial t} = - \frac{\partial J_w}{\partial l} - \frac{m a^{-1/b}}{l} u$$
This equation has the same form as the target equation if we identify $J_l = J_w$ and if the parameter $m$ in the target equation is understood as a new parameter $m_{new} = m_{old} a^{-1/b}$.

### **Step 3: Transform the Flux $J_w$**

Now, we express the flux $J_w$ in terms of $u$ and $l$. The expression for $J_w$ is:
$$J_w = (A w^n - B w) N - D w^{n+1} \frac{\partial N}{\partial w}$$
We need to substitute for $w$, $N$, and $\frac{\partial N}{\partial w}$.
First, let's find the derivative $\frac{\partial N}{\partial w}$ in terms of $u$ and $l$:
$$\frac{\partial N}{\partial w} = \frac{dl}{dw} \frac{\partial N}{\partial l} = \frac{1}{ab l^{b-1}} \frac{\partial}{\partial l} \left( \frac{u}{ab l^{b-1}} \right)$$
$$\frac{\partial N}{\partial w} = \frac{1}{ab l^{b-1}} \left[ \frac{1}{ab l^{b-1}}\frac{\partial u}{\partial l} - u \frac{(b-1)}{ab l^b} \right] = \frac{1}{(ab)^2 l^{2b-2}}\frac{\partial u}{\partial l} - \frac{u(b-1)}{(ab)^2 l^{2b-1}}$$
Now substitute this and the expressions for $w$ and $N$ into $J_w$:

* **First part of $J_w$ (advection term):**
    $$(A w^n - B w) N = (A (al^b)^n - B (al^b)) \frac{u}{ab l^{b-1}}$$
    Using $n=1-1/b \Leftrightarrow bn = b-1$:
    $$(A a^n l^{b-1} - B a l^b) \frac{u}{ab l^{b-1}} = \left(\frac{A a^{n-1}}{b} - \frac{B}{b} l\right) u$$

* **Second part of $J_w$ (diffusion term):**
    $$- D w^{n+1} \frac{\partial N}{\partial w} = - D (al^b)^{n+1} \left[ \frac{1}{(ab)^2 l^{2b-2}}\frac{\partial u}{\partial l} - \frac{u(b-1)}{(ab)^2 l^{2b-1}} \right]$$
    The exponent of $l$ is $b(n+1) = b(1-1/b+1) = 2b-1$.
    $$- D a^{n+1} l^{2b-1} \left[ \frac{1}{a^2b^2 l^{2b-2}}\frac{\partial u}{\partial l} - \frac{u(b-1)}{a^2b^2 l^{2b-1}} \right]$$
    $$= - \frac{D a^{n-1}}{b^2} l \frac{\partial u}{\partial l} + \frac{D a^{n-1}(b-1)}{b^2} u$$

Combining all parts of $J_w$:
$$J_w = \left(\frac{A a^{n-1}}{b} + \frac{D a^{n-1}(b-1)}{b^2}\right) u - \frac{B}{b} l u - \frac{D a^{n-1}}{b^2} l \frac{\partial u}{\partial l}$$

### **Step 4: Identify the Coefficients**

We now have the transformed PDE:
$$\frac{\partial u}{\partial t} = - \frac{\partial J_w}{\partial l} - \dots$$
with $J_w$ as derived above. We compare this flux with the target flux $J_l$:
$$J_l = k (L_\infty^0 - l) u - \alpha l \frac{\partial u}{\partial l} = (k L_\infty^0) u - k l u - \alpha l \frac{\partial u}{\partial l}$$
By setting $J_l = J_w$ and comparing the coefficients of the terms $u$, $lu$, and $l \frac{\partial u}{\partial l}$, we find the expressions for $k$, $L_\infty^0$, and $\alpha$.

* **Coefficient of $l \frac{\partial u}{\partial l}$:**
    $$-\alpha = - \frac{D a^{n-1}}{b^2} \Leftrightarrow \alpha = \frac{D a^{n-1}}{b^2}$$

* **Coefficient of $lu$:**
    $$-k = -\frac{B}{b} \Leftrightarrow k = \frac{B}{b}$$

* **Coefficient of $u$:**
    $$\begin{split}
    &k L_\infty^0 = \frac{A a^{n-1}}{b} + \frac{D a^{n-1}(b-1)}{b^2}\\
    \Leftrightarrow\qquad &L_\infty^0 = \frac{1}{k} \left(\frac{A a^{n-1}}{b} + \frac{D a^{n-1}(b-1)}{b^2}\right) = \frac{b}{B} a^{n-1} \left(\frac{A}{b} + \frac{D(b-1)}{b^2}\right)\\
    \Leftrightarrow\qquad &L_\infty^0 = \frac{a^{n-1}}{B} \left(A + \frac{D(b-1)}{b}\right)
    \end{split}$$

### **Final Expressions**

The expressions for $k$, $L_\infty^0$, and $\alpha$ in terms of $A$, $B$, $D$, and the transformation parameters $a$ and $b$ are:

$$\begin{split}
k &= \frac{B}{b},\\
\alpha &= \frac{D a^{n-1}}{b^2},\\
L_\infty^0 &= \frac{a^{n-1}}{B} \left(A + D\frac{b-1}{b}\right).
\end{split}$$

These can also be expressed in terms of $n$ by using the relations $n=1-1/b$, which implies $b = \frac{1}{1-n}$ and $n-1 = -1/b$. Also, the term $\frac{b-1}{b} = 1 - \frac{1}{b} = n$.

Substituting these into the expressions gives:

$$\begin{split}
k &= B(1-n),\\
\alpha &= D a^{-(1-n)} (1-n)^2,\\
L_\infty^0 &= \frac{a^{-(1-n)}}{B} (A + Dn).
\end{split}$$


**Note:** The parameter $L_\infty^0$ is related to the asymptotic length $L_\infty$ in
the expression for the average growth rate $$g(l) = k(L_\infty - l)$$ by
$$L_\infty = L_\infty^0 + \frac{\alpha}{k}.$$


### **Inverting the Relations**

We will now solve express the parameters $A$, $B$, and $D$ in terms of
$k$, $L_\infty^0$ and $\alpha$.

**1. Solving for B**

From the first equation, we can directly isolate $B$:
$$\begin{split}
&k = B(1-n)\\
\Leftrightarrow\qquad & B = \frac{k}{1-n}
\end{split}$$

**2. Solving for D**

From the second equation, we can directly isolate $D$:
$$\begin{split}
&\alpha = D a^{-(1-n)} (1-n)^2\\
\Leftrightarrow\qquad & D = \frac{\alpha}{a^{-(1-n)} (1-n)^2}\\
\Leftrightarrow \qquad & D = \frac{\alpha a^{1-n}}{(1-n)^2}
\end{split}$$

**3. Solving for A**

From the third equation, we can isolate $A$:
$$L_\infty^0 = \frac{a^{-(1-n)}}{B} (A + Dn)$$
Rearranging the terms to solve for $A$:
$$\begin{split}
&B L_\infty^0 = a^{-(1-n)} (A + Dn),\\
\Leftrightarrow \qquad &B L_\infty^0 a^{1-n} = A + Dn,\\
\Leftrightarrow \qquad &A = B L_\infty^0 a^{1-n} - Dn.
\end{split}$$
Now, we substitute the expressions we found for $B$ and $D$:
$$A = \left(\frac{k}{1-n}\right) L_\infty^0 a^{1-n} - \left(\frac{\alpha a^{1-n}}{(1-n)^2}\right) n$$
We can factor out common terms to simplify the expression for $A$:
$$A = \frac{a^{1-n}}{(1-n)} \left(k L_\infty^0 - \frac{n \alpha}{1-n}\right)$$

### **Final Inverted Relations**

The final expressions for $A$, $B$, and $D$ are:

$$\begin{split}
A &= \frac{a^{1-n}}{1-n} \left(k L_\infty^0 - \frac{n \alpha}{1-n}\right),\\
B &= \frac{k}{1-n},\\
D &= \frac{\alpha a^{1-n}}{(1-n)^2}.
\end{split}$$
