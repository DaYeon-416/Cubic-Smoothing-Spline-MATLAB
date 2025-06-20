# Cubic-Smoothing-Spline-MATLAB

## Files Used
- csaps_matlab.m
- Wage.csv (for example data)


## API
### Read Matrix
- xdata, ydata : 입력값으로 들어가는 데이터는 column vector로 이루어짐

    ```
    % For example...

    xdata = cell2mat(age_wage_array(:, 1)); % size : [61 x 1]
    ydata = cell2mat(age_wage_array(:, 2)); % size : [61 x 1]
    ```


### Parameter Selection
#### [weights]
- Values
    - ( 0, inf ) / NaN

        ```
        % For example...

        weights = ones(size(xdata)) * 0.1;
        weights = NaN; % NaN일 경우 prepare_data 함수를 통해 weight는 1로 이루어진 column vector가 됨
        ```

- Features 
    - 관찰된 실제 데이터와 cubic smoothing spline에 의해 예측된 데이터 사이의 잔차에 할당된 가중치를 의미함
    - 데이터를 잘 fitting 하는 것과 smooth curve를 갖는 것, 둘 사이의 trade-off를 조절함
    - 특정 구간에 대해 error measure weight가 크다면 smoothness를 희생하더라도 spline을 관찰 데이터와 비슷하게 fitting 하고자 하며, 반대로 error measure weight가 작다면 spline이 관찰 데이터에서 다소 벗어나더라도 smooth한 형태를 가지는 경향이 있음
    - 입력값으로 들어가는 데이터와 사이즈가 동일해야 함
    - 특정 구간의 weight의 변경 또한 가능
    
        ```
        % For example... 

        % 일부 구간의 weight 0.001, 나머지 구간 weight 1
        weights = ones(size(xdata)) * 1;
        weights(25:50) = 0.001;
        ```

---


#### [p]
- Values
    - [ 0, 1 ] / NaN

- Features
    - p가 NaN으로 지정되었다면 auto smoothing을 통해 p가 자동으로 계산됨
    - p가 0에 가까울수록 curve가 더욱 smooth 해지며, 1에 가까울수록 원본 데이터에 더욱 가깝게 fitting 됨

        ```
        % For example...

        p = 0;
        p = 0.001;
        p = 1;
        ```

- Parameter Selection 부분에서 어떤 파라미터를 사용하는지에 따라 최종적으로 결정되는 smoothing 파라미터 p는 변화함

    - if p = (0, 1), 
        - make_spline

    - if p = NaN, 
        - make_spline > compute_smooth

    - if normalizedsmooth = true, 
        - make_spline > normalize_smooth

---


#### [normalizedsmooth]
- Values
    - true / false

- Features
    - Function
        #### [normalize_smooth]

        - smoothing 파라미터를 통해 곡선의 smoothing 정도를 바꿀 수 있으나 이는 데이터의 규모에 따라 달라지기 때문에 데이터 정규화를 통해 이러한 문제를 해결할 수 있음
        - 데이터를 정규화한 후, 새로 정규화된 데이터를 기반으로 스무딩 파라미터가 계산되며 해당 매개변수를 사용하여 평활화 스플라인을 정규화하지 않고 원본 데이터에 맞출 수 있음
        - 즉, normalizedsmooth를 true로 변경 시, 정규화를 통해 smoothing 파라미터가 데이터의 영향을 받지 않고 데이터의 기본 패턴을 기반으로 파라미터가 선택될 수 있도록 만들어줌



### Function
#### [compute_smooth]
- Feature
    - smoothing 파라미터 p에 NaN값이 입력되었을 때 활성화됨
    - 식 (6(1-p) * QTD2Q + p * R) * c 에서 (6(1-p) * QTD2Q + p * R)은 coefficient를 이루는 행렬임
    - 계수의 행렬은 p * A + (1 - p) * B 의 형태를 띄고 있는데, p가 주어지지 않을 경우 p * A = (1 - p) * B이라고 가정함
    - 그러므로 p * trace(A) = (1 - p) * trace(B) 식을 풀면 p를 구할 수 있음

        ```math
        \begin{equation}
        \frac{p}{1 - p} = \frac {trace(B)}{trace(A)} = \frac {trace(6Q^T D^2 Q)}{trace(R)}
        \end{equation}
        ```
        ```math
        \begin{equation}
        \text{if} ~\frac{p}{1 - p} = r, ~~~p = \frac{r}{1+r}
        \end{equation}
        ```
---


#### [make_spline]
- Feature
    - x, y 데이터, weight, smoothing 파라미터, smoothing 정규화 여부에 따라 coefficient(계수)와 새로운 smoothing 파라미터 p를 계산함

```
dx = diff(x);
```
```math
h_i = x_{i+1} - x_i = \triangle{x_i}
```

```
diags_R  = vertcat(dx(2:end), 2*(dx(2:end) + dx(1:end-1)), dx(1:end-1));
R        = spdiags(diags_R', -1:1, N-2, N-2);
```
```math
\begin{equation*}
\underset{(n-2)\times(n-2)}R = 
\begin{bmatrix}
2(\triangle x_1+ \triangle x_2) &\triangle x_2 & 0 & \cdots & 0 & 0 \\
\triangle x_2 & 2(\triangle x_2+ \triangle x_3) & \triangle x_3 & \cdots &  0 & 0 \\
0 & \triangle x_3 & 2(\triangle x_3+ \triangle x_4)  & \cdots &  0 & 0 \\
\vdots  & \vdots  & \vdots  & \ddots & \vdots & \vdots\\
0 & 0 & 0 & \cdots & 2(\triangle x_{n-3}+ \triangle x_{n-2}) & \triangle x_{n-2}\\
0 & 0 & 0 & \cdots & \triangle x_{n-2} & 2(\triangle x_{n-2}+ \triangle x_{n-1})
\end{bmatrix}
\end{equation*}
```
    
```
diags_QT = vertcat(dx_denom(1:end-1), -(dx_denom(2:end) + dx_denom(1:end-1)), dx_denom(2:end));
QT       = spdiags(diags_QT', [0, 1, 2], N-2, N);
```
```math
\begin{equation*}
\underset{(n-2)\times(n)}{Q^T} = 
\begin{bmatrix}
\frac{1}{\triangle x_1} &  -(\frac{1}{\triangle x_1} + \frac{1}{\triangle x_2}) & \frac{1}{\triangle x_2} & 0 & \cdots  & 0& 0 & 0 \\
0 & \frac{1}{\triangle x_2} & -(\frac{1}{\triangle x_2} + \frac{1}{\triangle x_3})  & \frac{1}{\triangle x_3} & \cdots  & 0 & 0 & 0 \\
0 & 0 & \frac{1}{\triangle x_3}   & -(\frac{1}{\triangle x_3} + \frac{1}{\triangle x_4})& \cdots  & 0 & 0 & 0 \\
\vdots  & \vdots  & \vdots  & \ddots  &   & \vdots & \vdots & \vdots\\
0 & 0 & 0 & 0 & \cdots & -(\frac{1}{\triangle x_{n-3}} + \frac{1}{\triangle x_{n-2}}) & \frac{1}{\triangle x_{n-2}} & 0\\
0 & 0 & 0 & 0 & \cdots &  \frac{1}{\triangle x_{n-2}} & -(\frac{1}{\triangle x_{n-2}} + \frac{1}{\triangle x_{n-1}}) & \frac{1}{\triangle x_{n-1}}
\end{bmatrix}
\end{equation*}
```

```
D2 = spdiags((1./w)', 0, N, N);
```
```math
\begin{equation}
\underset{(n\times n)}D = 
\begin{bmatrix}
\delta y_1 &0 & 0 & \cdots & 0 \\
0 & \delta y_2 & 0 & \cdots &  0 \\
0 & 0 &  \delta y_3  & \cdots &  0  \\
\vdots  & \vdots  & \vdots  & \ddots & \vdots \\
0 & 0 & 0 & \cdots & \delta y_n  &\\
\end{bmatrix}
\end{equation}
```

```
a = 6.*(1.-p) * QTD2Q + p*R;
b = QT * y';
u = a \ b; % solve for u
```
```math
\begin{equation}
\underset{\mathbf{a}}{\left(6(1-p) Q^T D^2 Q+p R\right)} \mathbf{u}=\underset{\mathbf{b}}{Q^T \mathbf{y}}
\end{equation}
```

```
% cofficient 'a'

D2Qu = D2 * QT' * u;
yi   = y' - 6.*(1.-p) * D2Qu;
c4   = yi(1:end-1, :);
```
```math
\begin{equation}
\mathrm{a}_i=\mathrm{y}_i-6(1-p) D^2 Q \mathrm{u}_i
\end{equation}
```

```
% coefficient 'b'

c = 3 * p * u;
c_pad = zeros(size(c, 1) + 2, 1); % zero padding
c_pad(2:end-1) = c;
c = c_pad;

c3 = diff(yi)./dx - dx .* (2/3 * c(1:end-1, :) + 1/3 * c(2:end, :)); % b
```
```math
\begin{equation}
b_i=\left(a_{i+1}-a_i\right) \frac{1}{h_i}-\frac{2 c_i+c_{i+1}}{3} h_i
\end{equation}
```

```
% coefficient 'c'

c2 = c(1:end-1, :);
```
```math
\begin{equation}
\mathrm{c_i}=3 p \mathrm{u_i}
\end{equation}
```
```
% coefficient 'd'

c1 = diff(c) ./ (3 * dx); % d
```
```math
\begin{equation}
d_i = \frac{1}{3h_i}(c_{i+1}-c_i)
\end{equation}
```
---


#### [construct_fast]
```math
\begin{equation}
\text{for i = }1, 2, ..., {n-1},
\end{equation}
```
```math
\begin{equation}
S_{i}(x) = a_{i}+ b_{i}(x-x_{i})+ c_{i}(x-x_{i})^2+ d_{i}(x-x_{i})^3
\end{equation}
```
- Features
    - 관측 데이터인 break_point_list(csv의 x데이터), 절단점인 x_list, coefficient matrix를 통해 3차 다항식을 구성하기 위해 필요한 새로운 절단점 x에 대응하는 y값(y_list = S_i(x))을 구할 수 있음
    - nu
        - 0 : 미분 없음
        - 1 : 1차 미분식
        - 2 : 2차 미분식
