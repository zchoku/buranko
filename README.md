#buranko

運動方程式
$$
\begin{align}
r(t) = l - a\cos(\omega+\delta)\\
2 \dot{r}\dot{\theta} + r\ddot{theta} + g\sin(\theta)
\end{align}
$$
lはブランコの支点から中腰での重心の位置までの長さ、
aは立ち上がった状態での重心位置から中腰の状態での重心位置までの距離、
または中腰の状態での重心位置からしゃがんだ状態での重心位置までの距離
omegaは重心の運動の振動数
deltaは重心の運動の初期位相

codeでは
$$
\begin{equation}
\omega = k*\sqrt{\frac{g}{l}}
\end{equation}
$$
として、kを入力する。

また、t=0からt=6*Tまで、Nt分割で計算する。
Tは長さlの単振り子の場合の周期。プロット範囲を変更する場合はコードを書き換えて欲しい。

buranko_rungekutta_3
では、コンパイルして、
./a.out Nt delta k a
（Ntなどにはそれぞれ実数を入力して欲しい）
と実行すれば、与えたパラメータに応じた振幅の時間変化のグラフが表示される。
例えば、
./a.out 100 0. 2. 0.2
と実行すれば、分割数100、初期位相0. 振動数$$\omega = 2*\sqrt{g/l}$$,
重心の振幅a=0.2
での計算結果が表示される。


buranko_rungekutta_2.cpp
コンパイルして、
./a.out Nt delta k a
