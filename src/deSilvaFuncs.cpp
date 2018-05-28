#include <Rcpp.h>
using namespace Rcpp;

// Functions translated from the SAS code of de Silva et al (2005)
// doi: 10.1038/sj.hdy.6800728
// Used in deSilvaFreq and meandistance.matrix2.

// G function from de Silva et al.
// [[Rcpp::export]]
int G(int q, int n){
  int r = 1;
  for(int j = 0; j <= q; j++){
    r = r * (n + j) / (j + 1);
  }
  return r;
}

// INDEXG subroutine from de Silva et al.
// [[Rcpp::export]]
int INDEXG(IntegerVector ag1, int na1, int m2){
  int x = 1 + ag1[m2-1] - ag1[m2-2];
  for(int q = 1; q <= m2-2; q++){
    x = x + G(q, na1+1-ag1[m2-q-2]) - G(q, na1+1-ag1[m2-q-1]);
  }
  x = x + G(m2-1,na1) - G(m2-1,na1+1-ag1[0]);
  return x;
}

// GENLIST subroutine from de Silva et al.
// [[Rcpp::export]]
IntegerMatrix GENLIST(int ng, int na1, int m2){
  // set up temporary genotype vector and genotype array
  IntegerVector ag1(m2, 1);
  IntegerMatrix ag(ng, m2);
  
  // fill genotype array with all possible combinations
  int g = 0;
  for(int i = 0; i < m2; i++){
    ag(g, i) = ag1[i];
  }
  int a = m2 - 1;
  while(a > -1){
    if(ag1[a] == na1){
      a--;
    } else {
      if(a > -1){
        ag1[a]++;
        if(a < m2 - 1){
          for(int a1 = a+1; a1 < m2; a1++){
            ag1[a1] = ag1[a];
          }
        }
        g++;
        for(int i = 0; i < m2; i++){
          ag(g, i) = ag1[i];
        }
        a = m2 - 1;
        }
      }
    }
    return ag;
}

// RANMUL subroutine from de Silva et al.
// [[Rcpp::export]]
List RANMUL(int ng, int na1, IntegerMatrix ag, int m2){
  // rmul is multiplier to get genotype freq under random mating
  // arep shows how many copies of each allele each genotype has
  IntegerVector rmul(ng);
  IntegerMatrix arep(ng, na1);
  
  for(int g = 0; g < ng; g++){
    rmul[g] = 1;
    for(int i = 0; i < na1; i++){
      arep(g, i) = 0;
    }
    arep(g, ag(g,0) - 1) = 1;
    for(int j = 1; j < m2; j++){
      rmul[g] = rmul[g] * (j + 1);
      if(ag(g, j) == ag(g, j - 1)){
        arep(g, ag(g, j) - 1) = arep(g, ag(g, j) - 1) + 1;
        rmul[g] = rmul[g]/arep(g, ag(g, j) - 1);
      } else {
        arep(g, ag(g, j) - 1) = 1;
      }
    }
  }
  
  return List::create(Named("rmul") = rmul, Named("arep") = arep);
}

// SELFMAT subroutine from de Silva et al.
// [[Rcpp::export]]
IntegerMatrix SELFMAT(int ng, int na1, IntegerMatrix ag, int m2){
  // SELFMAT subroutine
  // smat is selfing matrix, and smatdiv is divisor for selfing matrix
  int m = m2/2;
  IntegerMatrix smat(ng, ng);
  IntegerVector al1(m);
  IntegerVector al2(m);
  IntegerVector ag1(m2);
  int a1;
  int a2;
  int j1;
  int j2;
  int k1;
  int k2;
  int g2;
  
  for(int g = 0; g < ng; g++){
    al1[0] = 1;
    if(m > 1){
      for(int j = 1; j < m; j++) al1[j] = al1[j - 1] + 1;
    }
    al1[m - 1]--;
    a1 = m - 1;
      while(a1 > -1){
        if(al1[a1] == (m + a1 + 1)){
          a1--;
        } else {
          if(a1 > -1){
            al1[a1]++;
            if(a1 < m - 1){
              for(int a3 = a1 + 1; a3 < m; a3++) al1[a3] = al1[a3 - 1] + 1;
            }
            al2[0] = 1;
            if(m > 1){
              for(int j = 1; j < m; j++) al2[j] = al2[j - 1] + 1;
            }
            al2[m - 1]--;
            a2 = m - 1;
            while(a2 > -1){
              if(al2[a2] == (m + a2 + 1)){
                a2--;
              } else {
                if(a2 > -1){
                  al2[a2]++;
                  if(a2 < m - 1){
                    for(int a3 = a2 + 1; a3 < m; a3++){
                      al2[a3] = al2[a3 - 1] + 1;
                    }
                  }
                  // UPDATESMAT subroutine
                  j1 = 0;
                  j2 = 0;
                  k1 = al1[j1];
                  k2 = al2[j2];
                  for(int i = 0; i < m2; i++){
                    if(k1 < k2){
                      ag1[i] = ag(g, k1 - 1);
                      if(j1 == m - 1){
                        k1 = 999;
                      } else {
                        j1++;
                        k1 = al1[j1];
                      }
                    } else {
                      ag1[i] = ag(g, k2 - 1);
                      if(j2 == m - 1){
                        k2 = 999;
                      } else {
                        j2++;
                        k2 = al2[j2];
                      }
                    }
                  }
                  g2 = INDEXG(ag1, na1, m2);
                  smat(g, g2 - 1)++;
                  // end UPDATESMAT subroutine
                  a2 = m - 1;
                }
              }
            }
            a1 = m - 1;
          }
        }
      }
  }
  return smat;
}
