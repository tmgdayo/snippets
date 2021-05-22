#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using vi = vector<int>;
using vl = vector<ll>;
using vvi = vector<vi>;
using vvl = vector<vl>;
#define REP(i, n) for (int i = 0, i##_len = (n); i < i##_len; ++i)
#define FOR(i, a, b) for (int i = (a), i##_len = (b); i < i##_len; i++)
#define ALL(x) (x).begin(), (x).end()  // sortなどの引数を省略したい
#define UNUSED_VARIABLE(x) ((void)(&x))  // 使わない引数の警告を消す
#ifdef _DEBUG
#define PRE_COMMAND             \
    std::cin.rdbuf(in.rdbuf()); \
    std::cout << fixed << setprecision(15);
#else
#define PRE_COMMAND                    \
    cout << fixed << setprecision(15); \
    ios::sync_with_stdio(false);       \
    cin.tie(0);
#endif
const double EPS = 1e-10, PI = acos(-1.0);
// 出力用関数
void PRINT() {}
template <class Head> void PRINT(Head&& head) { std::cout << head << '\n'; }
template <class Head, class... Tail> void PRINT(Head&& head, Tail&&... tail) {
    std::cout << head << ' ';
    PRINT(std::forward<Tail>(tail)...);
}
template <class T> auto MAX(T& seq) { return *max_element(ALL(seq)); }
template <class T> auto MIN(T& seq) { return *min_element(ALL(seq)); }
template <class T> auto SUM(T& seq) {
    T temp{0};
    auto& temp2 = temp[0];
    return accumulate(seq.begin(), seq.end(), temp2);
}
template <class T> void REV(vector<T>& seq) { reverse(ALL(seq)); }
template <class T> void SORT(T& seq) { sort(ALL(seq)); }
template <class T, class S> void SORT(T& seq, S& sort_order) {
    sort(ALL(seq), sort_order);
}
template <class T> void SORTR(vector<T>& seq) { sort(ALL(seq), greater<T>()); }
template <class T> int pow(int n_0, T k, int mod) {
    if (n_0 >= mod) n_0 %= mod;
    if (n_0 < 0) n_0 += mod;
    unsigned long long n = (unsigned long long)n_0, now = 1;
    while (k) {
        if (k & 1) now = (now * n) % mod;
        n = (n * n) % mod;
        k >>= 1;
    }
    return (int)now;
}
template <typename T> istream& operator>>(istream& is, vector<T>& vec) {
    for (T& x : vec) is >> x;
    return is;
}
template <typename T> ostream& operator<<(ostream& os, const vector<T>& v) {
    if (!v.size()) return os;
    typename vector<T>::const_iterator ii = v.begin();
    os << *ii++;
    for (; ii != v.end(); ++ii) os << " " << *ii;
    return os;
}
void yn(bool flag) {
    if (flag) {
        PRINT("YES");
    } else {
        PRINT("NO");
    }
}