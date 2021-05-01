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
bool prime[1000001];

void Quadratic_residue(int n);
void Legendre_symbol();
int Legendre_symbol(int a, int p);
void primes_sieve();
map<int, int> prime_factors(int n);
void Jacobi_symbol();

signed main() {
	fast;
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
#endif

	int n = 11;
	// cin<<n;
	// Quadratic_residue(n);
	// Legendre_symbol();
	Jacobi_symbol();
	// cout << Legendre_symbol(-35, 97);

	// primes_sieve();
	// for (int p : primes) {
	// 	if (Legendre_symbol(10, p) == 1)
	// 		cout << "FOUNDDDDD: " << p;
	// }

	return 0;
}

void Quadratic_residue(int p) {
	cout << "## QUADRATIC RESIDUES of " << p << "\n\n";
	// input is a prime
	set<int> QR, QNR;
	for (int i = 1; i < p; i++) {
		QR.insert(power(i, 2, p));
	}
	for (int i = 1; i < p; i++) {
		if (QR.find(i) == QR.end())
			QNR.insert(i);
	}
	cout << "Quadratic residues:\n";
	for (int a : QR) {
		cout << a << " ";
	}
	cout << "\nQuadratic NON residues:\n";
	for (int a : QNR) {
		cout << a << " ";
	}
}

void Legendre_symbol() {
	int a = -23, p = 83;
	// cin>>a>>p;
	a = (a + p) % p;
	cout << "## Legendre Symbol of " << a << "/" << p << "\n\n";
	int L;

	if (a == 0) {
		L = 0;
	}
	if (a == p - 1) {
		if ((p % 4) == 1)
			L = 1;
		else L = -1;
	}
	else if (a == 2) {
		if ((p % 8) == 1 || (p % 8) == 7)
			L = 1;
		else L = -1;
	}

	else {
		L = power(a, (p - 1) / 2, p);
		if (L == p - 1)
			L = -1;
	}

	cout << " Legendre_symbol value is: " << L << "\n\n";
}

int Legendre_symbol(int a, int p) {
	cout << "\n## Legendre Symbol of " << a << "/" << p << "\n";
	int L;
	a = (a + p) % p;
	cout << " = (" << a << "/" << p << ")\n\n";

	if (a == 0) {
		L = 0;
	}
	if (a == p - 1) {
		if ((p % 4) == 1)
			L = 1;
		else L = -1;
	}
	else if (a == 2) {
		if ((p % 8) == 1 || (p % 8) == 7)
			L = 1;
		else L = -1;
	}

	else {
		cout << " power step: " << a << "^" << (p - 1) / 2 << " mod" << p << "\n";
		L = power(a, (p - 1) / 2, p);
		if (L == p - 1)
			L = -1;
	}

	cout << " Legendre_symbol value is: " << L << "\n";
	return L;
}

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

map<int, int> prime_factors(int n) {
	primes_sieve();
	// finding the prime factors of n
	int k = 0;
	map<int, int> factors;
	for (auto i : primes) {
		if (i > n)	break;
		if (n % i == 0 && n > 1) {
			k = 0;
			while (n % i == 0 && n > 1) {
				n /= i;
				k++;
			}
			factors[i] = k;
		}
	}
	return factors;
}

void Jacobi_symbol() {
	int a = 11, n = 61;
	cin >> a >> n;
	a = (a + n) % n;

	if (gcd(a, n) != 1) {
		cout << " GCD of (a,n) != 1, implies Jacobi symbol = 0 \n";
		return;
	}

	cout << "## Jacobi Symbol of " << a << "/" << n << "\n\n";
	int J = 1;
	// factorize n
	cout << "factors of n : \nn = ";
	auto mp = prime_factors(n);
	for (auto t : mp) {
		cout << t.first << "^" << t.second << " * ";
	}
	for (auto t : mp) {
		J *= (power(Legendre_symbol(a, t.first), t.second));
	}
	cout << "\n Jacobi_symbol value is: " << J << "\n\n";

}