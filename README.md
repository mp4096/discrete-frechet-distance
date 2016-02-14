## Discrete Frechét distance
Compute the discrete Frechét distance [1] between two curves specified by ordered discrete points in _n_-dimensional space according to algorithm in [2].

The implementation is based on MATLAB function by Zachary Danziger [3].
This function is 10 to 50 times as fast as [3] at the cost of computing only the DFD (no coupling sequence).

### Installation
Open the repo directory in MATLAB and type
```matlab
mex DiscreteFrechetDistance.c
```
If you want to add the directory to your path, type
```matlab
addpath(cd)
savepath
```

### Usage example
:warning: In contrast to [3], this function requires that input matrices contain point coordinates in _columns_, i.e. the matrix shape should be `<number of dimensions> x <number of points in curve {1, 2}>`.

```matlab
x1 = linspace(-2, 2, 200);
y1 = x1.^2 + 0.2.*x1 + 2;

t3 = linspace(0, 4, 300);
x2 = 8.*sin(0.1.*t3) - 1;
y2 = 2.*cos(0.9.*t3 + 0.2) + 6;

t3 = linspace(0, 1, 1000);
x3 = 2.*sin(10.*t3);
y3 = 10.*cos(2.*t3 + 0.2) + 3;

figure
plot(x1, y1);
hold on
plot(x2, y2);
plot(x3, y3);
hold off
legend('C1', 'C2', 'C3')

% Notice the matrix layout
d12 = DiscreteFrechetDistance([x1; y1], [x2; y2]);
d23 = DiscreteFrechetDistance([x3; y3], [x2; y2]);
d13 = DiscreteFrechetDistance([x3; y3], [x1; y1]);

fprintf('DFD(C1, C2) = %5.2f\n', d12);
fprintf('DFD(C2, C3) = %5.2f\n', d23);
fprintf('DFD(C3, C1) = %5.2f\n', d13);
```

### Known issues
* Microsoft Windows SDK 7.1 seems to have problems with `fmin` and `fmax` from `math.h`. Compile this function using MinGW 4.9.2 C/C++ (TDM-GCC).

### References:
[1] Wikipedia: [Frechét distance](https://en.wikipedia.org/wiki/Fr%C3%A9chet_distance)

[2] Eiter, T. and Mannila, H.: Computing Discrete Frechét Distance. Technical Report 94/64, TU Wien, 1994; available [\[online\]](http://www.kr.tuwien.ac.at/staff/eiter/et-archive/cdtr9464.pdf).

[3] Danziger, Z.: Discrete Frechét Distance. Available [\[online\]](http://www.mathworks.com/matlabcentral/fileexchange/31922-discrete-frechet-distance).
