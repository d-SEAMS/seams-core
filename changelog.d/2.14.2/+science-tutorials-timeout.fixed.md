The science tutorial check keeps a 180 s timeout. An AddressSanitizer build runs that check serially instead of allowing 900 s: under contention the tutorial took 212 s, and alone it took 22.6 s.
