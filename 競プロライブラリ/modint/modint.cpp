template <unsigned int mod> struct modint {
    unsigned int val;
    // _mod < 2^31である必要あり
    // _modが素数でないときは割り算ができない
    modint() : val(0) {}
    template <class T> modint(T v) {
        if (v >= 0) {
            if (v < (T)mod) {
                val = (unsigned int)(v);
            } else {
                val = (unsigned int)(v % (long long)mod);
            }
        } else {
            if (-v > (T)mod) {
                val = (unsigned int)(v + mod);
            } else {
                val = (unsigned int)(v % (long long)mod + mod);
            }
        }
    }
    // 単項演算子++の定義
    modint& operator++() {
        val++;
        if (val == mod) val = 0;
        return *this;
    }
    // 単項演算子--の定義
    modint& operator--() {
        if (!val) val = mod;
        val--;
        return *this;
    }
    // 二項演算子++の定義
    modint operator++(int) {
        modint result = *this;
        ++*this;
        return result;
    }
    // 二項演算子--の定義
    modint operator--(int) {
        modint result = *this;
        --*this;
        return result;
    }
    // +=の定義
    modint& operator+=(const modint& rhs) {
        val += rhs.val;
        if (val >= mod) val -= mod;
        return *this;
    }
    // -=の定義
    modint& operator-=(const modint& rhs) {
        if (val < rhs.val) val += mod;
        val -= rhs.val;
        return *this;
    }
    // *=の定義
    modint& operator*=(const modint& rhs) {
        unsigned long long v = val;
        v *= rhs.val;
        val = (unsigned int)(v % mod);
        return *this;
    }
    // /=の定義
    modint& operator/=(const modint& rhs) { return *this = *this * rhs.inv(); }
    // 単項演算子boolの定義
    explicit operator bool() const { return val; }
    // 単項演算子-の定義
    bool operator!() const { return !val; }
    // 二項演算子+の定義
    friend modint operator+(const modint& lhs, const modint& rhs) {
        return modint(lhs) += rhs;
    }
    // 二項演算子-の定義
    friend modint operator-(const modint& lhs, const modint& rhs) {
        return modint(lhs) -= rhs;
    }
    // 二項演算子*の定義
    friend modint operator*(const modint& lhs, const modint& rhs) {
        return modint(lhs) *= rhs;
    }
    // 二項演算子/の定義
    friend modint operator/(const modint& lhs, const modint& rhs) {
        return modint(lhs) /= rhs;
    }
    // 二項演算子==の定義
    friend bool operator==(const modint& lhs, const modint& rhs) {
        return lhs.val == rhs.val;
    }
    // 二項演算子!=の定義
    friend bool operator!=(const modint& lhs, const modint& rhs) {
        return lhs.val != rhs.val;
    }
    // 単項演算子+の定義
    modint operator+() const { return *this; }
    // 単項演算子-の定義
    modint operator-() const { return modint() - *this; }
    // cinによる入力
    friend istream& operator>>(istream& is, modint& rhs) {
        long long v;
        istream& ret = is >> v;
        if (v >= 0) {
            if (v < mod) {
                rhs.val = (unsigned int)(v);
            } else {
                rhs.val = (unsigned int)(v % (long long)mod);
            }
        } else {
            if (-v > mod) {
                rhs.val = (unsigned int)(v + mod);
            } else {
                rhs.val = (unsigned int)(v % (long long)mod + mod);
            }
        }
        return ret;
    }
    // coutによる出力
    friend ostream& operator<<(ostream& os, const modint& rhs) {
        return os << rhs.val;
    }
    // powの定義、使う演算子は既にmodint内で定義しているのでintのときと同様でOK
    modint pow(long long n) const {
        modint x(val), r(1);
        while (n) {
            if (n & 1) r *= x;
            x *= x;
            n >>= 1;
        }
        return r;
    }
    // 逆元を出力
    modint inv() const {
        modint x(val), r(1);
        unsigned int n = mod - 2;
        while (n) {
            if (n & 1) r *= x;
            x *= x;
            n >>= 1;
        }
        return r;
    }
    // modを出力
    unsigned int mod_() const { return mod; }
};