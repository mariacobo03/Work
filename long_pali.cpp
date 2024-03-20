# include <iostream>
# include <map>
# include <string>
# include <set>

using namespace std;


bool palindrome(const string& s) {
	map <char, char> d = {{'A','T'}, {'C','G'}, {'G','C'}, {'T','A'}};
	int low = 0;
	int high = s.length()-1;
	while(low < high){
		if (d[s[low]] != s[high]){
			return false;
		}
		else {
			low = low + 1;
			high = low - 1;
		}
	}
	return low != high;
}

set<string> longest_palindromic(const string& s){
	set<string> seq;
	int big = 0;
	for (int i=0; i < int(s.length()) - 1; ++ i){
		for (int j = i + 1; j < int(s.length()); ++ j){
			if (j - i + 1 >= big and palindrome(s.substr(i, j - i + 1))){
				if (j - i + 1 == big){
					seq.insert(s.substr(i, j - i + 1));
				}
				else {
					seq.clear();
					seq.insert(s.substr(i, j - i + 1));
					big = j - i + 1;
				}
			}
		}
	}
	return seq;
}


int main(){
    string s;
    cin >> s;
    set<string> seq = longest_palindromic(s);
    for(auto x : seq) cout << x << endl;
}

