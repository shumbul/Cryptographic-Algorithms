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

vi primes;
bool prime[100001];

void primes_sieve() {
	// cout << "PRIME SEIVE\n";
	memset(prime, true, sizeof prime);
	int n = sizeof(prime) - 1;
	for (int p = 2; p * p <= n ; p++) {
		if (prime[p] == true) {
			for (int i = p * p; i <= n; i += p)
				prime[i] = false;
		}
	}
	for (int p = 2; p <= n; p++)
		if (prime[p])
			primes.pb(p);

	// for0(i, primes.size()) {
	// 	cout << primes[i] << " -> ";
	// }
}

int phi(int n) {
	primes_sieve();
	int res = n;
	if (prime[n]) {return n - 1;}
	for (int i = 2; i * i <= n; i++) {
		if (n % i == 0) {
			res /= i, res *= (i - 1);
			while (n % i == 0) n /= i;
		}
	}
	if (n > 1)	res /= n, res *= (n - 1);
	return res;
}


bool iscong(int x, int b, int m) {
	return ((b - x) % m == 0);
}

void Linear_cong() {
	cout << " LINEAR CONGRUENCE \n ";
	int t;	cin >> t;
	while (t--) {
		int a, b, m, M, A, B;
		// ax = bmodm
		cin >> a >> b >> m;
		M = m;
		A = a, B = b;
		// has a soluion only when gcd(a,m)|b
		int g = gcd(a, m);
		// cout << g;
		if (b % g) {
			cout << "GCD of a, m does not divide v, NO SOLUTION!\n";
			return;
		}
		a /= g, b /= g, m /= g;
		cout << "ax cong b mod m = > " << a << "x cong " << b << " mod " << m << "\n";
		//x = (a^(phi(m)-1))*b mod m
		int phi_m = phi(m);
		cout << "\n phi " << phi_m << "\n";
		int apow = power(a, (phi_m - 1), m);
		int x = (apow) * b % m;
		cout << "x = " << a << "^" << phi_m - 1 << " * " << b << " mod " << m << "\n";
		cout << "x = " << apow << " * " << b << " mod " << m << "\n";
		cout << "x = " << apow * b << " mod " << m << "\n";

		cout << "\nSolutions: \n";
		while (x <= M) {
			cout << x << ", ";
			if (iscong((A * x), B, M) == 0)
				cout << "Test FAILED!!!\n";
			x += m;
		}
		if (iscong((A * x), B, M)) {
			cout << "Test passed\n\n";
		}
		else cout << "Test FAILED!!!\n\n";
	}
}

signed main() {
	fast;
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
#endif

	Linear_cong();	//t, i=1 to t: {a, b, m}	ax cong bmodm

	return 0;
}