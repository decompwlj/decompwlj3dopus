\\ decompwlj.gp — the reference kernels, and a CSV writer for the 3-D atlas.
\\
\\ Project convention: the reference implementation is the readable one.  The
\\ C generator in this folder computes the same function far faster; these are
\\ what it is checked against.
\\
\\ Verified on PARI/GP 2.15.4.
\\   gp -q tools/decompwlj.gp

\\ ---------------------------------------------------------------- kernels

\\ The definition, transcribed.  O(tau(l)) after one factorisation.
decomp_fact(n, n1) = {
  my(d, l, D);
  d = n1 - n; l = n - d;
  if (l <= d, return([0, 0, d]));          \\ not decomposable
  D = divisors(l);                         \\ sorted ascending; D[#D] = l > d
  for (i = 1, #D, if (D[i] > d, return([D[i], l/D[i], d])));
}

\\ The compact fordiv form (decompwlj_fordiv.txt).
decomp(a, b) = {
  my(d = b - a, l);
  if (a <= 2*d, return([0, 0, d]));
  l = a - d;
  fordiv(l, k, if (k > d, return([k, l/k, d])));
}

\\ Classification of a result triple: 1 level (k > L), 0 weight (k <= L,
\\ ties included), -1 unclassified (not decomposable).
dclass(t) = if (t[1] == 0, -1, if (t[1] > t[2], 1, 0));

\\ ------------------------------------------------------- the atlas format
\\
\\ One chunk file per CHUNK rows, header "a,d,k,L", EVERY term written —
\\ a non-decomposable one as k = L = 0 — so that row i of the concatenation
\\ is always sequence index n = n0 + i.  (Census engines that skip
\\ non-decomposable terms silently shift every index that is keyed on n.)

CHUNK = 50000;

\\ v: a vector of terms, strictly increasing, length >= 2.  The last term has
\\ no successor in v and so is not written; ask for one term more than you want.
write_chunks(dir, v) = {
  my(N = #v - 1, c = 0, f, t);
  system(Str("mkdir -p ", dir));
  while (c * CHUNK < N,
    f = Str(dir, "/chunk-", if(c<10,"00",if(c<100,"0","")), c, ".csv");
    write(f, "a,d,k,L");
    for (i = c*CHUNK + 1, min((c+1)*CHUNK, N),
      t = decomp_fact(v[i], v[i+1]);
      write(f, v[i], ",", t[3], ",", t[1], ",", t[2]);
    );
    c++;
  );
  printf("%s: %d rows in %d chunk(s)\n", dir, N, c);
}

\\ Example — the first 2000 primes into data/seq/primes_gp/
\\   write_chunks("data/seq/primes_gp", vector(2001, i, prime(i)))

\\ ----------------------------------------------------------------- checks
\\ Table 1 of the treatise, and the deep witnesses.
{
  my(ok = 1);
  forprime (p = 2, 10^5,
    if (decomp_fact(p, nextprime(p+1)) != decomp(p, nextprime(p+1)), ok = 0));
  printf("kernels agree below 10^5: %d\n", ok);
  print(decomp_fact(1000000000039, nextprime(1000000000040)));  \\ [461, 2169197397, 22]
  print(decomp_fact(10^18 + 3, nextprime(10^18 + 4)));          \\ [47, 21276595744680851, 6]
  print(decomp_fact(10^24 + 7, nextprime(10^24 + 8)));          \\ [11909, 83970106642035435385, 42]
}
