// importing libraries

# include <set>
# include <string>
# include <iostream>
using namespace std;

string gene_finding(string s);

int main() {
    string s;
    cin >> s;
    cout << gene_finding(s) << endl;
}

// afterwards we define the gene_finding function, in pyhton does not matter but in C++ yes

string gene_finding(string s) {
    set<string> stop = {"TAA", "TAG", "TGA"};
    int n = s.length();
    int i = 0;
    while ( i + 2 <= n) {
        if (s.substr(i, 3) == "ATG") {
            int j = i + 3;
            while(j + 2 <= n) {
                if (stop.find(s.substr(j, 3)) != stop.end())
                    return s.substr(i, j + 2 - i + 1);
                j = j + 3;
        } }
        i = i + 1;
    }
    return "";
}

/* s.substr(j, 3) --> substring of position j of size 3 
stop.end() --> it goes one more positon when it stops
if it is different from end it means it didn't find it */




