# Tesi triennale Binelli

Questo programma fa parte della tesi triennale di Riccardo Binelli dal 
titolo "Il decadimento alpha come fenomeno di tunneling quantistico:
modello a potenziale e simulazioni numeriche", depositato presso la
Biblioteca delle Scienze dell'Università degli Studi di Pavia 
nell'Anno Accademico 2024/2025.
Il programma è stato sviluppato per calcolare i tempi di dimezzamento
degli isotopi pari dell'Uranio (dal 224U al 238U) sfruttando un modello
a potenziale di tipo Coulomb + Woods-Saxon e una routine di integrazione
della famiglia delle regole di Gauss-Kronrod, di cui viene qui presentato
il codice e le cui caratteristiche si trovano all'interno del testo di
tesi.

All'interno del file `main.cc` sono presenti i parametri da noi utilizzati,
così che una volta compilato l'eseguibile restituisca la tabella dei tempi
simulati. Nel caso si volesse richiamare la funzione di integrazione 
`integrate` all'interno di altri progetti, lo usage è il seguente:
```c pp
integrate(a, b, function, accuracy, depth)
```
dove [a,b] è l'intervallo di integrazione e dove accuracy e depth sono 
settati di default se non specificati. 

I risultati sono riportati e discussi all'interno della tesi.
