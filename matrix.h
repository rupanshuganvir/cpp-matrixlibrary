#include<iostream>
using namespace std;
template<class T>
class matrix{
// m: Number of rows, n: NUmber of columns
    public:
    int m;
    int n;
    T **element;

    matrix(int r, int c){
        m=r; n=c;
        element = new T*[r];
        for (int i=0; i<r; i++) {element[i]=new T[c];}
    }
    matrix(int r, int c, T **input){
        m=r; n=c;
        element = new T*[r];
        for (int i=0; i<r; i++) {
            element[i]=new T[c];
            for (int j=0; j<c; j++){ element[i][j] = input[i][j];}
        }
    }
    matrix(int r, int c, T input[]){
        m=r; n=c;
        int k=0;
        element = new T*[r];
        for (int i=0; i<r; i++) {
            element[i]=new T[c];
            for (int j=0; j<c; j++){ element[i][j] = input[k]; k++; }
        }
    }
    void swap (int r1, int r2){
        T *temp =element[r1];
        element[r1]=element[r2];
        element[r2]=temp;
        return;
    }
    // 1 -> lamda*r
    void op1 (int r, T lamda){
        for (int i=0; i<n; i++) element[r][i] *=lamda;
        return;
    }
    // r1 -> r1 + lamda*r2
    void op2 (int r1, int r2, T lamda){
        for (int i=0; i<n; i++) element[r1][i] += lamda*element[r2][i];
        return;
    }
    matrix<T> operator=(const matrix<T> &M){
        return matrix<T>(M.m,M.n,M.element);
    }
    matrix<T> operator+(const matrix<T> &M){
        matrix<T> ans=M;
        for (int i=0; i<M.m; i++){
            for (int j=0; j <M.n; j++){
                ans.element[i][j]+= element[i][j];
            }
        }
        return ans;
    }
    matrix<T> operator-(const matrix<T> &M){
        matrix<T> ans=M;
        for (int i=0; i<M.m; i++){
            for (int j=0; j <M.n; j++){
                ans.element[i][j]-= element[i][j];
            }
        }
        return ans;
    }
    //scalar multiplication
    matrix<T> operator*(const T &lamda){
        return lamda*(*this);
    }
};

//the overloaded operator <<
template <class T>
ostream &operator<<(ostream &ost, const matrix<T> &X){
    for (int i=0; i<X.m; i++){
        ost<<"( ";
        for (int j=0; j<X.n; j++){ ost<<X.element[i][j]<<" "; }
        ost<<")\n";
    }
    return ost;
}

//to check if two matrices are equal.
//ISSUE: Check for template types as well.
template<class T>
bool operator==(const matrix<T> &A, const matrix<T> &B){
    if (A.m != B.m || A.n != B.n ) return false;
    for (int i=0; i<A.m; i++){
        for (int j=0; j<A.n; j++){
            if (A.element[i][j] != B.element[i][j]) return false;
        }
    }
    return true;
}

//Assuming that the matrices are multiplication compatible
template<class T>
matrix<T> operator*(const matrix<T> &A, const matrix<T> &B){
    matrix<T> ans(A.m, B.n);
    for (int i=0; i<A.m; i++){
        for (int j=0; j<B.n; j++){
            ans.element[i][j] =0;
            for (int k=0; k<A.n; k++) { ans.element[i][j] += A.element[i][k] * B.element[k][j]; }
        }
    }
    return ans;
}

//~A is the transpose of A
template<class T>
matrix<T> operator ~(const matrix<T> &M){
    matrix<T> ans(M.n,M.m);
    for (int i=0; i<M.n; i++){
        for (int j=0; j<M.m; j++){
            ans.element[i][j]=M.element[j][i];
        }
    }
    return ans;
}

//reduced form of a matrix
template<class T>
matrix<T> reduced(const matrix<T> &A){
    matrix<T> ans=A;
    // i: index of pivot row, j: index of pivot column
    for (int i=0,j=0; i<ans.m && j<ans.n; i++){
        if (ans.element[i][j]==0){
            int q=ans.m-1;
            while (ans.element[q][j]==0 && q!=i){q--;}
            if (q==i) { i--; j++; continue;}
            T* temp = ans.element[i];
            ans.element[i] = ans.element[q];
            ans.element[q]=temp;
        }
        else{
            for (int p=i+1; p<ans.m; p++){
                T factor = ans.element[p][j];
                for (int q=j; q<ans.n; q++){
                    ans.element[p][q] *= ans.element[i][j];
                    ans.element[p][q] -= ans.element[i][q]*factor;
                }
            }
            j++;
        }
    }
    for (int i=ans.m-1; i>0; i--){
        int j=0; while (ans.element[i][j]==0 && j<ans.n) {j++;} if (j==ans.n) {continue;}

        T lamda=ans.element[i][j]; for (int q=j; q<ans.n; q++) {ans.element[i][q] /= lamda;}
        for (int p=i-1; p>=0; p--) {
            lamda = ans.element[p][j];
            for (int q=j; q<ans.n; q++) { ans.element[p][q] -= ans.element[i][q]*lamda; }
        }
    }
    return ans;
}

template<class T>
matrix<T> C_space(const matrix<T> &A){
    matrix<T> R = reduced(~A);
    int r=0; //r:rank of A
    for (int i=0,j=0; i<R.m && j<R.n; ){
        if (R.element[i][j]!=0) {r++; i++; j++;}
        else{ j++; }
    }
    T *e = new T[r*R.n];
    int k=0;
    for (int i=0; i<r; i++){
        for (int j=0; j<R.n; j++){
            e[k++]=R.element[i][j];
        }
    }
    return ~matrix<T>(r, R.n, e);
}
