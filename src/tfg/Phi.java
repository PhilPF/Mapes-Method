/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package tfg;

import java.util.*;

/**
 *
 * @author ppita
 */
public class Phi {
    
    private final static double boundC = Math.log(12/Math.E);
    private int a;
    private long ma2;
    private long phima=1;
    private long[] table;
    
    public Phi(int a){
        this.a=a;
        generateTable();
    }
    
    private long[] table(){
        return table;
    }
    
    public long valueAtX(long x){
        if (x==0) return 0;
        if (a==0) return x;
        long absx=Math.abs(x);
        long r=absx%(2*ma2);
        long z=(absx/(2*ma2))*phima;
        //System.out.println(z);
        long c=table.length;
        if (r<ma2){
            for (int i=0; i<table.length; i++) {
                if((r+1)<=table[i]){
                    c=i;
                    break;
                }
            }
            z+=c;
        }  
        else {
            for (int i=0; i<table.length; i++) {
                if((2*ma2-r)<=table[i]){
                    c=i;
                    break;
                }
            }
            z+=phima-c;
        }
        return (long)(absx/x)*z;
    }
    
    private void generateTable(){
        long[] firstAPrimes=firstNPrimes(a);
        ma2=1;
        for (int i=1; i<a; i++) ma2*=firstAPrimes[i];
        if(a==0 || a==1) table=new long[]{1};
        else{
            table = TFG.eratosthenesInterval((int)firstAPrimes[firstAPrimes.length-1],ma2, firstAPrimes);
            table[0]=1;
        }
        for (long p:firstAPrimes) phima*=(p-1);
    }
    
    private long[] firstNPrimes(int n){
        
        int stop=(int) (12*(n*(Math.log(n)+boundC)));
        //System.out.println(stop);

        boolean[] isPrime = new boolean[stop];
        for(int i=0; i<stop; i++) isPrime[i]=true;

        int primeCount=1;
        for (int k=1; k<stop && primeCount<n; k++){
            if (isPrime[k]){
                int p=2*k+1, p2=p*p;
                if(p2<=stop){
                    for(int i=(p2-1)/2; i<stop; i+=p){
                        isPrime[i]=false;
                    }
                }
                primeCount++;
            }
        }
        
        int count=1;
        long[] newPrimes =new long[primeCount];
        newPrimes[0]=2;
        for(int i=1; i<stop && count<primeCount; i++){
            if (isPrime[i]) newPrimes[count++]=(2*i+1);
        } 
        
        return newPrimes;

    }
    
    
}
