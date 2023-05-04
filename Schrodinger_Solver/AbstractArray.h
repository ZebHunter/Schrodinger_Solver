#ifndef ARRAYT_HPP  

#define ARRAYT_HPP  



#include <cstdlib>
#include <cstring>
#include <iostream> 

using namespace std;



template < class T >
class arrayt {
public:
    
    arrayt(const int n1 = 1);              
    arrayt(const int n1, const int n2);   
    arrayt(const arrayt<T>& a);

    
    inline ~arrayt() { if (nn > 0) delete[] p; nn = nndim = nn1 = nn2 = 0; }

    
    inline arrayt<T>& operator=(const arrayt<T>& m);
    inline arrayt<T>& operator=(const T);
    inline T& operator()(const int i1, const int i2);      
    inline T& operator()(const int i);                    
    inline arrayt<T>& operator+=(const arrayt<T>& m);

   
    inline int n1() const { return nn1; }
    inline int n2() const { return nn2; }
    inline int ndim() const { return nndim; }
    inline int n() const { return nn; }
    void resize(const int n1, const int n2);   
    void resize(const int n);                  

private:    

    T* p;                   
    int nn;                 
    int nndim;             
    int nn1, nn2;           
};



template < class T >
arrayt<T>::arrayt(const int n1)           
{
    if (n1 <= 0) {
        cout << "arrayt initialized with size = " << n1 << ", NOT ALLOWED" << endl;
        exit(EXIT_FAILURE);
    }
    p = new T[n1];  
    nn = nn1 = n1;
    nn2 = 0;
    nndim = 1;
}

template < class T >
arrayt<T>::arrayt(const int n1, const int n2)     
{
    if ((n1 <= 0) || (n2 <= 0)) {
        cout << "arrayt initialized with size = "
            << n1 << " x " << n2 << ", NOT ALLOWED" << endl;
        exit(EXIT_FAILURE);
    }
    p = new T[n1 * n2];  
    nn1 = n1;
    nn2 = n2;
    nn = n1 * n2;
    nndim = 2;
}

template < class T >                
arrayt<T>::arrayt(const arrayt<T>& a)
{
    p = new T[a.nn];   
    memcpy(p, a.p, a.nn * sizeof(T));
    nn1 = a.nn1;
    nn2 = a.nn2;
    nndim = a.nndim;
    nn = a.nn;
}



template < class T >
void arrayt<T>::resize(const int n)    
{
    if (n <= 0) {
        cout << "arrayt resize() with size = " << n
            << ", NOT ALLOWED" << endl;
        exit(EXIT_FAILURE);
    }
    if (nn > 0) delete[] p;
    p = new T[n];   
    nn1 = n;
    nn2 = 0;
    nndim = 1;
    nn = n;
}

template < class T >
void arrayt<T>::resize(const int n1, const int n2)       
{
    if ((n1 <= 0) || (n2 <= 0)) {
        cout << "arrayt resize() with size = " << n1 << " x "
            << n2 << ", NOT ALLOWED" << endl;
        exit(EXIT_FAILURE);
    }
    if (nn > 0) delete[] p;
    p = new T[n1 * n2];  
    nn1 = n1;
    nn2 = n2;
    nndim = 2;
    nn = n1 * n2;
}





template < class T >
arrayt<T>& arrayt<T>::operator=(const arrayt<T>& m)
{
    if ((nn != m.nn) || (nndim != m.nndim)) {
        cout << "arrayt = operator invoked with unequal sizes\n"
            << "  m1 size = " << nn << ", dim= " << nndim << "\n"
            << "  m2 size = " << m.nn << ", dim= " << m.nndim << endl;
        exit(EXIT_FAILURE);
    }
    else {
        memcpy(p, m.p, nn * sizeof(T)); // fastest way to do this
        return *this;
    }
}


template < class T >
arrayt<T>& arrayt<T>::operator=(const T value)
{
    if (nn > 0)
        for (int i = 0; i < nn; i++) p[i] = value;
    return *this;
}


template < class T >
inline T& arrayt<T>::operator()(const int i1, const int i2)
{
#ifdef ARRAYT_BOUNDS_CHECK
    if ((i1 < 0) || (i1 >= nn1) ||
        (i2 < 0) || (i2 >= nn2) || (nndim != 2)) {
        cout << "out of bounds index in arrayt\n"
            << "  size = " << nn1 << " x " << nn2 << ", ndim= " << nndim << "\n"
            << "  access = (" << i1 << ", " << i2 << ")" << endl;
        exit(EXIT_FAILURE);
    }
#endif

 
    return *(p + i2 + i1 * nn2);       

}


template < class T >
inline T& arrayt<T>::operator()(const int i)
{
#ifdef ARRAYT_BOUNDS_CHECK
    if ((i < 0) || (i >= nn) || (nndim != 1)) {
        cout << "out of bounds index in arrayt\n"
            << "  size = " << nn << ", ndim= " << nndim << "\n"
            << "  access = " << i << endl;
        exit(EXIT_FAILURE);
    }
#endif

    return *(p + i);
}



template < class T >
inline arrayt<T>& arrayt<T>::operator+=(const arrayt<T>& m)
{
    if ((m.nn != nn)) {
        cout << "arrayt += operator invoked with unequal sizes:\n"
            "   " << nn << " and " << m.n() << endl;
        exit(EXIT_FAILURE);
    }
    else {
        int i;
        for (i = 0; i < nn; i++) p[i] += m.p[i];
        return *this;
    }
}

#endif  // ARRAYT_HPP