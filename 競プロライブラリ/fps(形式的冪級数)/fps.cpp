// 最大次数
template <size_t mxm, class mint> struct fps {
    using sfps = vector<tuple<mint, size_t>>;
    size_t n, SIZE_T_MAX = (size_t)(-1);  // 現在の次数
    unsigned int mod = mint().mod_();
    vector<mint> v;
    // コンストラクタ
    fps(vector<mint> tmp) {
        v = tmp;
        n = v.size() - 1;
        resize();
    }
    fps(mint tmp) {
        v = {tmp};
        n = 0;
    }
    fps() {
        v = {0};
        n = 0;
    }
    // resize、szは最大サイズ
    void resize(size_t sz = mxm) {
        if (sz < n) {
            n = sz;
            v.resize(n + 1);
        }
        if (n) {
            size_t i = n;
            for (; i and !v[i]; i--) {}
            n = i;
            v.resize(i + 1);
        }
    }
    // sparse用resize
    void sp_resize(const sfps& rhs) {
        if (n < get<1>(rhs.back())) {
            n = get<1>(rhs.back());
            v.resize(n + 1);
        }
    }
    // z = aを代入したときの値
    mint eval(const mint& a) {
        mint ret = 0, now = 1;
        for (mint& val : v) {
            ret += val * now;
            now *= a;
        }
        return ret;
    }
    // z^dの係数
    mint coef(const size_t& d) {
        if (d > n) return mint(0);
        return v[d];
    }
    // (1+c*(z^d))を乗じる
    fps& multiply(const mint& c, const size_t& d) {
        n = min(mxm, n + d);
        v.resize(n + 1);
        for (size_t i = n - d; i != SIZE_T_MAX; i--) v[i + d] += c * v[i];
        return *this;
    }
    // (1+c*(z^d))で割る
    fps& divide(const mint& c, const size_t& d) {
        n = mxm;
        v.resize(n + 1);
        for (size_t i = d; i <= n; i++) v[i] -= c * v[i - d];
        return *this;
    }
    // NTTを使った乗算
    // 逆元を計算する、szは最大サイズ
    fps inv(size_t sz = mxm) const {
        assert(v[0]);
        if (sz < 230) {
            fps f({(mint)1});
            vector<mint> fv(sz + 1, 0);
            fv[0] = 1;
            mint iv = (v[0]).inv();
            for (size_t i = 0; i <= sz; i++) {
                size_t d = min(n, i);
                mint tmp = fv[i];
                for (size_t j = 1; j <= d; j++) tmp -= fv[i - j] * v[j];
                fv[i] = tmp * iv;
            }
            return fps(move(fv));
        }
        size_t m = 1;
        vector<mint> res = {v[0].inv()};
        while (m <= sz) {
            vector<mint> f(v.begin(), v.begin() + min(n + 1, 2 * m));
            vector<mint> g = res;
            f.resize(2 * m), g.resize(2 * m);
            dft(f), dft(g);
            for (size_t i = 0; i < 2 * m; i++) f[i] *= g[i];
            idft(f);
            f.erase(f.begin(), f.begin() + m);
            f.resize(2 * m);
            dft(f);
            for (size_t i = 0; i < 2 * m; i++) f[i] *= g[i];
            idft(f);
            mint iz = mint(2 * m).inv();
            iz *= -iz;
            for (size_t i = 0; i < m; i++) f[i] *= iz;
            res.insert(res.end(), f.begin(), f.begin() + m);
            m *= 2;
        }
        res.resize(sz + 1);
        return fps(move(res));
    }
    // 微分を計算する
    fps diff() const {
        if (n == 0) return fps();
        vector<mint> f(n);
        for (size_t i = 1; i <= n; i++) f[i - 1] = i * v[i];
        size_t i = n - 1;
        for (; i and !f[i]; i--) {}
        f.resize(i + 1);
        return fps(move(f));
    }
    // 積分を計算する
    fps integral() const {
        vector<mint> v1(n + 2, 1), v2(n + 2, 0), f(n + 2);
        for (size_t i = 1; i <= n; i++) v1[i + 1] = v1[i] * i;
        v2[n + 1] = (v1[n + 1] * (n + 1)).inv();
        for (size_t i = n + 1; i > 1; i--) v2[i - 1] = v2[i] * i;
        for (size_t i = 0; i < n + 2; i++) v1[i] *= v2[i];
        for (size_t i = 1; i <= n + 1; i++) f[i] = v[i - 1] * v1[i];
        size_t i = n + 1;
        for (; i and !f[i]; i--) {}
        f.resize(i + 1);
        return fps(move(f));
    }
    // logを計算する
    fps log(size_t sz = mxm) const {
        assert(v[0] == 1);
        fps f = diff() * inv(sz);
        f.resize(sz);
        f = f.integral();
        f.resize(sz);
        return f;
    }
    // expを計算する
    fps exp(size_t sz = mxm) const {
        if (!n) return fps(1);
        vector<mint> g{1}, g_fft{1, 1}, f = v;
        f[0] = 1;
        f.resize(sz + 1);
        // 微分の前計算
        vector<mint> h_drv(v.begin() + 1, v.end());
        for (size_t i = 0; i < n; i++) h_drv[i] *= (i + 1);
        h_drv.resize(sz);
        // 積分用の前計算
        vector<mint> v1(2 * sz + 2, 1), v2(2 * sz + 2, 0);
        for (size_t i = 1; i <= 2 * sz; i++) v1[i + 1] = v1[i] * i;
        v2[2 * sz + 1] = (v1[2 * sz + 1] * (2 * sz + 1)).inv();
        for (size_t i = 2 * sz + 1; i > 1; i--) v2[i - 1] = v2[i] * i;
        for (size_t i = 0; i < 2 * sz + 2; i++) v1[i] *= v2[i];
        // 本体
        for (size_t m = 2; m <= sz; m *= 2) {
            // prepare
            vector<mint> f_fft(f.begin(), f.begin() + m);
            f_fft.resize(2 * m), dft(f_fft);
            mint m2 = mint(2 * m).inv();
            mint m1 = m2 * 2;
            mint m3 = -m1;
            mint m4 = m3 * m1;
            // Step 2.a'
            {
                vector<mint> _g(m);
                for (size_t i = 0; i < m; i++) _g[i] = f_fft[i] * g_fft[i];
                idft(_g);
                _g.erase(_g.begin(), _g.begin() + m / 2);
                _g.resize(m), dft(_g);
                for (size_t i = 0; i < m; i++) _g[i] *= g_fft[i];
                idft(_g);
                _g.resize(m / 2);
                for (size_t i = 0; i < m / 2; i++) _g[i] *= m4;
                g.insert(g.end(), _g.begin(), _g.begin() + m / 2);
            }
            // Step 2.b'--d'
            vector<mint> t(f.begin() + 1, f.begin() + m);
            for (size_t i = 0; i < m - 1; i++) t[i] *= (i + 1);
            t.resize(m);
            {
                // Step 2.b'
                vector<mint> r(h_drv.begin(), h_drv.begin() + m - 1);
                // Step 2.c'
                r.resize(m);
                dft(r);
                for (size_t i = 0; i < m; i++) r[i] *= f_fft[i];
                idft(r);
                for (size_t i = 0; i < m; i++) r[i] *= m3;
                // Step 2.d'
                for (size_t i = 0; i < m; i++) t[i] += r[i];
                t.insert(t.begin(), t.back());
                t.pop_back();
            }
            // Step 2.e'
            if (2 * m <= sz) {
                t.resize(2 * m);
                dft(t);
                g_fft = g;
                g_fft.resize(2 * m);
                dft(g_fft);
                for (size_t i = 0; i < 2 * m; i++) t[i] *= g_fft[i];
                idft(t);
                t.resize(m);
                for (size_t i = 0; i < m; i++) t[i] *= m2;
            } else {
                vector<mint> g1(g.begin() + m / 2, g.end());
                vector<mint> s1(t.begin() + m / 2, t.end());
                t.resize(m / 2);
                t.resize(m), dft(t);
                g1.resize(m), dft(g1);
                s1.resize(m), dft(s1);
                for (size_t i = 0; i < m; i++)
                    s1[i] = g_fft[i] * s1[i] + g1[i] * t[i];
                idft(s1);
                for (size_t i = 0; i < m; i++) t[i] *= g_fft[i];
                idft(t);
                for (size_t i = 0; i < m / 2; i++) t[i + m / 2] += s1[i];
                for (size_t i = 0; i < m; i++) t[i] *= m1;
            }
            // Step 2.f'
            vector<mint> u(f.begin() + m, f.begin() + min(sz + 1, 2 * m));
            t.insert(t.begin(), m - 1, 0);
            t.insert(t.begin(), 0);
            for (size_t i = 1; i < 2 * m; i++) t[i] *= v1[i];
            for (size_t i = 0; i < min(sz - m + 1, m); i++) u[i] -= t[m + i];
            // Step 2.g'
            u.resize(2 * m);
            dft(u);
            for (size_t i = 0; i < 2 * m; i++) u[i] *= f_fft[i];
            idft(u);
            for (size_t i = 0; i < m; i++) u[i] *= m2;
            // Step 2.h'
            for (size_t i = 0; i < min(sz + 1 - m, m); i++) f[m + i] = u[i];
        }
        return fps(move(f));
    }
    // powを計算する
    fps pow(size_t m, size_t sz = mxm) const {
        if (n == 0 and v[0] == 0) return fps(0);
        size_t l = 0;
        for (size_t i = 0; i < n; i++) {
            if (v[i] == 0) {
                l++;
            } else {
                break;
            }
        }
        fps g(v);
        g >>= l;
        if (g.n == 0 and g.v[0] == 0) return fps(0);
        mint c_inv = (g.v[0]).inv(), c_pow = (g.v[0]).pow(m);
        g *= c_inv;
        g = g.log(sz);
        g *= m;
        g = g.exp(sz);
        g *= c_pow;
        g <<= m * l;
        g.resize();
        return g;
    }
    // NTTを使った乗算
    fps& fast_multiply(const fps& rhs) {
        size_t exp = 1, sz = 2;
        for (; sz < n + rhs.n + 1; sz <<= 1, exp++) {}
        vector<mint> g(rhs.v);
        v.resize(sz), g.resize(sz);
        dft(v), dft(g);
        for (size_t i = 0; i < sz; i++) v[i] *= g[i];
        idft(v);
        mint iz = mint(sz).inv();
        for (size_t i = 0; i < sz; i++) v[i] *= iz;
        n = v.size() - 1;
        resize();
        return *this;
    }
    // NTTを使った除算
    fps& fast_divide(const fps& rhs) {
        fast_multiply(rhs.inv());
        n = v.size() - 1;
        resize();
        return *this;
    }
    // +=の定義
    fps& operator+=(const fps& rhs) {
        if (n < rhs.n) {
            n = rhs.n;
            v.resize(n + 1);
        }
        for (size_t i = 0; i <= rhs.n; i++) v[i] += rhs.v[i];
        resize();
        return *this;
    }
    // +=の定義（sparse）
    fps& operator+=(const sfps& rhs) {
        sp_resize(rhs);
        for (const tuple<mint, size_t>& t : rhs) v[get<1>(t)] += get<0>(t);
        resize();
        return *this;
    }
    // +=の定義（スカラー）
    fps& operator+=(const mint& rhs) {
        v[0] += rhs;
        return *this;
    }
    // -=の定義
    fps& operator-=(const fps& rhs) {
        if (n < rhs.n) {
            n = rhs.n;
            v.resize(n + 1);
        }
        for (size_t i = 0; i <= rhs.n; i++) v[i] -= rhs.v[i];
        resize();
        return *this;
    }
    // -=の定義（sparse）
    fps& operator-=(const sfps& rhs) {
        sp_resize(rhs);
        for (const tuple<mint, size_t>& t : rhs) v[get<1>(t)] -= get<0>(t);
        resize();
        return *this;
    }
    // -=の定義（スカラー）
    fps& operator-=(const mint& rhs) {
        v[0] -= rhs;
        return *this;
    }
    // *=の定義
    fps& operator*=(const fps& rhs) {
        size_t lhs_sz = n, rhs_sz = rhs.n;
        // 配列が大きい場合はfastで
        if (min(lhs_sz, rhs_sz) > 50) return fast_multiply(rhs);
        // 片方の配列が小さい場合はnaiveで
        n = min(mxm, lhs_sz + rhs_sz);
        v.resize(n + 1);
        for (size_t i = n; i != SIZE_T_MAX; i--) {
            mint tmp = 0;
            size_t d = 0, d2 = min(i, rhs.n);
            if (i > lhs_sz) d = i - lhs_sz;
            for (size_t j = d; j <= d2; j++) tmp += v[i - j] * rhs.v[j];
            v[i] = tmp;
        }
        resize();
        return *this;
    }
    // *=の定義（sparse）
    fps& operator*=(const sfps& rhs) {
        size_t lhs_sz = n, rhs_sz = get<1>(rhs.back());
        n = min(mxm, lhs_sz + rhs_sz);
        v.resize(n + 1);
        size_t k = rhs.size() - 1, d;
        mint m;
        for (size_t i = n; i != SIZE_T_MAX; i--) {
            if (k < 0) break;
            mint tmp = 0;
            for (size_t j = k; j != SIZE_T_MAX; j--) {
                tie(m, d) = rhs[j];
                if (i < d) {
                    k--;
                    continue;
                }
                if (d + lhs_sz < i) break;
                tmp += v[i - d] * m;
            }
            v[i] = tmp;
        }
        resize();
        return *this;
    }
    // *=の定義（スカラー）
    fps& operator*=(const mint& rhs) {
        for (size_t i = 0; i <= n; i++) v[i] *= rhs;
        resize();
        return *this;
    }
    // /=の定義
    fps& operator/=(const fps& rhs) {
        size_t rhs_sz = rhs.n;
        // 分母の配列が大きい場合はfastで
        if (min(n, rhs_sz) > 300) return fast_divide(rhs);
        // 分母の配列が小さい場合はnaiveで
        n = mxm;
        v.resize(n + 1);
        assert((rhs.v)[0]);
        mint rhs_inv = ((rhs.v)[0]).inv();
        for (size_t i = 0; i <= n; i++) {
            size_t d = min(rhs_sz, i);
            mint tmp = v[i];
            for (size_t j = 1; j <= d; j++) tmp -= v[i - j] * rhs.v[j];
            v[i] = tmp * rhs_inv;
        }
        resize();
        return *this;
    }
    // /=の定義（sparse）
    fps& operator/=(const sfps& rhs) {
        size_t rhs_sz = rhs.size();
        n = mxm;
        v.resize(n + 1);
        assert(get<0>(rhs[0]) and !get<1>(rhs[0]));
        mint rhs_inv = get<0>(rhs[0]).inv(), m;
        size_t d;
        for (size_t i = 0; i <= n; i++) {
            mint tmp = v[i];
            for (size_t j = 1; j < rhs_sz; j++) {
                tie(m, d) = rhs[j];
                if (i < d) break;
                tmp -= v[i - d] * m;
                if (j == rhs_sz - 1) break;
            }
            v[i] = tmp * rhs_inv;
        }
        resize();
        return *this;
    }
    // /=の定義（スカラー）
    fps& operator/=(const mint& rhs) {
        for (size_t i = 0; i <= n; i++) v[i] /= rhs;
        resize();
        return *this;
    }
    // >>=の定義
    fps& operator>>=(const size_t& rhs) {
        if (rhs == 0) return *this;
        if (rhs > n) {
            n = 0;
            v = {0};
        } else {
            for (size_t i = rhs; i <= n; i++) v[i - rhs] = v[i];
            n -= rhs;
            v.resize(n + 1);
        }
        return *this;
    }
    // <<=の定義
    fps& operator<<=(const size_t& rhs) {
        if (rhs == 0) return *this;
        if (rhs > mxm) {
            n = 0;
            v = {0};
            return *this;
        }
        n = min(n + rhs, mxm);
        v.resize(n + 1);
        for (size_t i = n; i >= rhs; i--) v[i] = v[i - rhs];
        for (size_t i = 0; i < rhs; i++) v[i] = 0;
        resize();
        return *this;
    }
    // =の定義
    fps operator=(vector<mint>&& rhs) {
        n = rhs.size() - 1;
        this->v = rhs;
        return *this;
    }
    // 二項演算子+の定義
    template <class T> friend fps operator+(const fps& lhs, const T& rhs) {
        return fps(lhs) += rhs;
    }
    // 二項演算子-の定義
    template <class T> friend fps operator-(const fps& lhs, const T& rhs) {
        return fps(lhs) -= rhs;
    }
    // 二項演算子*の定義
    template <class T> friend fps operator*(const fps& lhs, const T& rhs) {
        return fps(lhs) *= rhs;
    }
    // 二項演算子/の定義
    template <class T> friend fps operator/(const fps& lhs, const T& rhs) {
        return fps(lhs) /= rhs;
    }
    // 二項演算子>>の定義
    friend fps operator>>(const fps& lhs, const size_t& rhs) {
        return fps(lhs) >>= rhs;
    }
    // 二項演算子<<の定義
    friend fps operator<<(const fps& lhs, const size_t& rhs) {
        return fps(lhs) <<= rhs;
    }
    // 二項演算子==の定義
    friend bool operator==(const fps& lhs, const fps& rhs) {
        return lhs.v == rhs.v;
    }
    // 二項演算子!=の定義
    friend bool operator!=(const fps& lhs, const fps& rhs) {
        return lhs.v != rhs.v;
    }
    // 単項演算子+の定義
    fps operator+() const { return *this; }
    // 単項演算子-の定義
    fps operator-() const {
        fps ret = v;
        for (mint& m : ret.v) m = -m;
        return ret;
    }
    // 高速フーリエ変換を行う
    void dft(vector<mint>& x) const {
        size_t xs = x.size();
        size_t exp = 1, tmp = 2;
        for (; tmp < xs; tmp <<= 1, exp++) {}
        mint e = mint(3).pow((mod - 1) >> exp), now = 1;
        size_t m = xs >> 1;
        vector<mint> exp_list(m);
        for (mint& el : exp_list) {
            el = now;
            now *= e;
        }
        for (size_t cnt = 0; cnt < exp; cnt++) {
            for (size_t i = 0; i < xs; i += 2 * m) {
                mint* p = exp_list.data();
                for (size_t j = i; j < i + m; j++) {
                    mint u = x[j], t = x[j + m];
                    x[j] = u + t;
                    x[j + m] = (u - t) * (*p);
                    p++;
                }
            }
            m >>= 1;
            for (size_t i = 0; i < m; i++) exp_list[i] = exp_list[2 * i];
        }
    }
    // 高速逆フーリエ変換を行う
    void idft(vector<mint>& x) const {
        size_t xs = x.size(), exp = 1, tmp = 2;
        for (; tmp < xs; tmp <<= 1, exp++) {}
        mint e = mint(3).pow((mod - 1) >> exp).inv();
        vector<mint> exp_list(exp);
        for (mint& el : exp_list) {
            el = e;
            e *= e;
        }
        size_t m = 1;
        for (size_t cnt = exp - 1;; cnt--) {
            mint now = 1, w = exp_list[cnt];
            vector<mint> exp_list2(m);
            for (mint& el : exp_list2) {
                el = now;
                now *= w;
            }
            for (size_t i = 0; i < xs; i += 2 * m) {
                mint* p = exp_list2.data();
                for (size_t j = i; j < i + m; j++) {
                    mint u = x[j], t = x[j + m] * (*p);
                    x[j] = u + t;
                    x[j + m] = u - t;
                    p++;
                }
            }
            m <<= 1;
            if (!cnt) break;
        }
    }
    // cinによる入力
    friend istream& operator>>(istream& is, fps& rhs) { return is >> rhs.v; }
    // coutによる出力
    friend ostream& operator<<(ostream& os, const fps& rhs) {
        return os << rhs.v;
    }
};