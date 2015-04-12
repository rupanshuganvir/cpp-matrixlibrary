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
    matrix(int r, int c, T lamda){
        m=r, n=c;
        element = new T*[r];
        for (int i=0; i<r; i++) {
            element[i]=new T[c];
            for (int j=0; j<c; j++){ if (i==j) element[i][j]=lamda; else element[i][j]=0;}
        }
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
int rank(const matrix<T> &A){
    matrix<T> R = reduced(A);
    int r=0; //r:rank of A
    for (int i=0,j=0; i<R.m && j<R.n; ){
        if (R.element[i][j]!=0) {r++; i++; j++;}
        else{ j++; }
    }
    return r;
}

template<class T>
matrix<T> C_space(const matrix<T> &A){
    matrix<T> R = reduced(~A);
    int r=rank(R);
    T *e = new T[r*R.n];
    int k=0;
    for (int i=0; i<r; i++){
        for (int j=0; j<R.n; j++){
            e[k++]=R.element[i][j];
        }
    }
    return ~matrix<T>(r, R.n, e);
}

template<class T>
matrix<T> diag(T *diagonal, int n){
    T **e = new T*[n];
    for (int i=0; i<n; i++) e[i] = new T[n];
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            if (i==j) e[i][j] = diagonal[i];
            else e[i][j]=0;
        }
    }
    return matrix<T>(n,n,e);
}

//Asumming A to be square
template<class T>
T det(const matrix<T> &A){
    matrix<T> R=A;
    T det=1;
    for (int i=0; i<R.n;i++){
        if (R.element[i][i]==0) {
            int q=i+1;
            while (R.element[q][i]==0 && q<R.n) {q++;}
            if (q==R.n) return 0;
            R.swap(i,q); det*=-1;
            i--; continue;
        }
        else{
            for (int p=i+1; p<R.n; p++){
                T lamda = - (R.element[p][i]/R.element[i][i]);
                R.op2(p,i, lamda );
            }
        }
    }
    for (int i=0; i<R.n; i++){ det *= R.element[i][i]; }
    return det;
}

template<class T>
matrix<T> inverse(const matrix<T> &A){
    matrix<T> X=A; matrix<T> ans(X.n,X.n,1);
    for (int i=0; i<X.n;i++){
        if (X.element[i][i]==0) {
            int q=i+1;
            while (X.element[q][i]==0 && q<X.n) {q++;}
            if (q==X.n) {/*this will imply that its noninvertible*/ }
            X.swap(i,q); ans.swap(i,q);
            i--; continue;
        }
        else{
            for (int p=i+1; p<X.n; p++){
                T lamda = 1/X.element[i][i];
                X.op1(i,lamda); ans.op1(i,lamda);
                lamda = -(X.element[p][i]/X.element[i][i]);
                X.op2(p,i,lamda ); ans.op2(p,i,lamda);
            }
        }
    }
    for (int i=A.n-1; i>0;i--){
        for (int p=i-1; p>=0; p--){
            T lamda = -X.element[p][i];
            X.op2(p,i,lamda); ans.op2(p,i,lamda);
        }
    }
    return ans;
}

