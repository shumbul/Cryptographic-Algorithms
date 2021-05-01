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

int inv(int a, int m) {
	int m0 = m, t, q;
	int x0 = 0, x1 = 1;
	if (m == 1)
		return 0;
	// Apply extended Euclid Algorithm
	int num = 99999;
	while (a > 1 && num >= 0) {
		// q is quotient
		q = a / m;
		t = m;
		// m is remainder now, process same as
		// euclid's algo
		m = a % m, a = t;
		t = x0;
		x0 = x1 - q * x0;
		x1 = t;
		num--;
	}
	if (num == 0) {
		cout << "Infinite while loop\n";
		return 0;
	}
	// Make x1 positive
	if (x1 < 0)
		x1 += m0;
	return x1;
}

void findMinX(int num[], int rem[], int k) {
	// Compute product of all numbers
	int prod = 1;
	for (int i = 0; i < k; i++)
		prod *= num[i];
	cout << "Product of num's = " << prod << "\n\n";
	cout << "mi\tbi\t\txi\t\tNi\t\tbixiNi\n";
	// Initialize result
	int result = 0;
	// Apply above formula
	for (int i = 0; i < k; i++) {
		int pp = prod / num[i];
		int xi = inv(pp, num[i]);
		int xbN = rem[i] * xi * pp;
		result += xbN;
		cout << num[i] << "\t" << rem[i] << "\t\t" << xi << "\t\t" << pp << "\t\t" << xbN << "\n";
	}
	cout << "Summation (xibiNi) = " << result << "\n";
	int x = result % prod;
	cout << "\nx is: " << x << " mod " << prod << "\n\n";

	for0(i, k) {
		if (iscong(x, rem[i], num[i]) == 0)
			cout << "TEST FAILED at " << i << " !!\n";
		cout << "Passed: " << x << " " << rem[i] << " " << num[i] << "\n";
	}
}

void Chinese_rem() {
	cout << "CHINESE REMAINDER THEOREM\n\n";
	int n;	cin >> n;	//number of equations
	int num[n], rem[n];
	for0(i, n)	cin >> rem[i] >> num[i];
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if (gcd(num[i], num[j]) != 1) {
				cout << "Chinese remainder theorem not valid in this case!\n As each par of modulo should be relatively prime for it to hold.\n";
				cout << num[i] << " and " << num[j] << " are not co - primes\n";
				return;
			}
		}
	}
	findMinX(num, rem, n);

	for0(i, n) {
		cout << "Phi " << num[i] << " -> " << phi(num[i]) << "\n";
	}
	/*
	Nixi = 1 mod ni
	xi = Ni^(phi(ni)-1) mod ni
	*/
}

void Chinese_rem(int n, int rem[], int num[]) {
	cout << "CHINESE REMAINDER THEOREM after conversions\n\n";
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if (gcd(num[i], num[j]) != 1) {
				cout << "Chinese remainder theorem not valid in this case!\n As each par of modulo should be relatively prime for it to hold.\n";
				cout << num[i] << " and " << num[j] << " are not co - primes\n";
				return;
			}
		}
	}
	findMinX(num, rem, n);

	for0(i, n) {
		cout << "Phi " << num[i] << " -> " << phi(num[i]) << "\n";
	}
	/*
	Nixi = 1 mod ni
	xi = Ni^(phi(ni)-1) mod ni
	*/
}

void Chinese_converted() {
	cout << "Equations conversions - CHINESE REMAINDER\n\n";
	int n;	cin >> n;
	int a[n], b[n], m[n];
	for0(i, n)	cin >> a[i] >> b[i] >> m[i];
	int new_b;
	cout << "b's and m's for chinese remainder\n";
	cout << n << "\n";
	for0(i, n) {
		new_b = (power(a[i], phi(m[i]) - 1, m[i]) * b[i]) % m[i];
		b[i] = new_b;
		cout << b[i] << " " << m[i] << "\n";
	}
	Chinese_rem(n, b, m);
}

signed main() {
	fast;
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
#endif

	Chinese_rem();	//n, i=1 to n: {b[i], m[i]}
	// Chinese_converted();	//n, i=1 to n: {a[i], b[i], m[i]}


	return 0;
}