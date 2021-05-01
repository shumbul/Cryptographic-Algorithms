#include<bits/stdc++.h>
using namespace std;

typedef long long ll;
#define int long long
#define vi vector<int>
#define pb push_back
#define mp make_pair
#define ff first
#define ss second
#define inar(arr,n) for(int i=0; i<n; i++) cin>>arr[i];
#define suma(arr,n) for(int i=0; i<n; i++) sum+=arr[i];
#define for0(i,n) for(int i=0; i<n; i++)
#define for1(i,n) for(int i=1; i<=n; i++)
#define forr(i,n) for(int i=n-1; i>=0; i--)
#define rof(i,a,b) for(int i=a;i<=b;i++)
#define all(v) v.begin(), v.end()
#define rall(v) v.rbegin(), v.rend()
#define print(v) cout<<v<<"\n"
#define fast ios_base::sync_with_stdio(false); cin.tie(0); cout.tie(0)
const int mod = 1e9 + 7;
ll gcd (ll a, ll b) {return ( a ? gcd(b % a, a) : b );}
ll power(ll a, ll n) {ll p = 1; while (n > 0) {if (n % 2) {p = p * a;} n >>= 1; a *= a;} return p;}
ll power(ll a, ll n, ll mod) {ll p = 1; while (n > 0) {if (n % 2) {p = p * a; p %= mod;} n >>= 1; a *= a; a %= mod;} return p % mod;}
int n;
void Pollards_rho();
void Pollards_pminus1();

signed main() {
	fast;
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
#endif

	n = 387058387;
	cin >> n;

	Pollards_rho();
	// Pollards_pminus1();

	return 0;
}

int f(int x) {
	return (power(x, 2) + 1 ) % n;
}

void Pollards_rho() {
	cout << " ### POLLARDS RHO Method: for n = " << n << "\n";
	cout << "Step-1: choose a non-linear map function \n(fill above function)\n";

	cout << "\nStep-2: perform successive iteration k times \n(choose x0 and k below)\n";
	int x0 = 3, k = 4;
	map<int, int> X;
	X[0] = x0;
	cout << "x0 = " << x0 << "\n";
	for (int i = 1; i <= k; i++) {
		X[i] = f(X[i - 1]);
		cout << "x" << i << " = " << X[i] << "\n";
	}

	cout << "\nStep-3: compute g = gcd(xj-xk,n)\n";
	cout << "(if g is a non-trivial factor value, it's also a non-trivial factor of n\n";
	for (int i = 0; i < k; i++) {
		for (int j = i + 1; j <= k; j++) {
			int g = gcd(abs(X[i] - X[j]), n);
			cout << " GCD(abs(x" << i << " - x" << j << "), " << n << ") = " << g << "\n";
			if (g != 1 && g != n) {
				cout << "\n ### We got a non-trivial factor!!\n";
				cout << g << " is a non trivial factor of " << n << "\n";
				cout << " Another factor is " << n / g << "\n\n";
				return;
			}
		}
	}
	cout << "\n Could not find a non-trivial factor for n\n";
}

int lcm(int a, int b)
{
	return (a / gcd(a, b)) * b;
}

int lcm_till(int B) {
	int ans = 2;
	for (int i = 3; i <= B; i++) {
		ans = lcm(ans, i);
	}
	return ans;
}

void Pollards_pminus1() {
	// for d > 1 if m has a prime factor p such that (p - 1) divides (n!)
	cout << " ### POLLARDS p-1 Method: for n = " << n << "\n";
	cout << "Step-1: choose integer B>2 (fill below)\n";
	int B = 9;
	cout << " B = " << B << "\n";

	cout << "\n Step-2: choose k, which is a multiple of most of the numbers b<=B\n";
	int k = lcm_till(B);
	cout << " k = " << k << "\n";

	cout << "\n Step-3: choose a, which is a random number 2<=a<=n-2\n";
	int a = 2 + rand() % (n - 2);
	a = 2;	// for simplicity!!
	cout << " a = " << a << "\n";

	cout << "\n Step-4: compute r, r cong to a^k modn\n";
	int r = power(a, k, n);
	cout << " r = " << r << "\n";

	cout << "\n Step-5: compute d, d=gcd(r-1, n)\n";
	int d = gcd(r - 1, n);
	cout << " d = " << d << "\n";

	if (d == 1 || d == n)
		cout << " => Go to Step-1\n";
	else {
		cout << "\n ### We got a non-trivial factor!!\n";
		cout << d << " is a non trivial factor of " << n << "\n";
		cout << " Another factor is " << n / d << "\n\n";
	}
}