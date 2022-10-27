#include <iostream>
#include <vector>

using namespace std;

int main() {
	vector<int> seq = {0,1,4,9,16,25,36,49,64,81,100};
	for (const int &i : seq) {
		cout << i << " ";
	}
	cout << endl;

	for (int i = 1; i < (int) seq.size(); i++) {
		seq[i] += seq[i-1];
	}	

	for (const int &i : seq) {
		cout << i << " ";
	}
	cout << endl;
	return 0;
}
