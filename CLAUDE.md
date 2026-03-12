I want to test the comprehensiveness of Wycheproof's test vectors. Mutation testing is a nice way to do it, but for something like a signature verifier there's a risk you'll skip a check entirely in the implementation and then there won't be anything to mutate to notice test vectors for it are missing. 

If one can procure an implementation where the bugs are likely to be uncorrelated from one's own, then mutation testing that against your vectors can figure out gaps. 

Make an ML-DSA verifier (no need to implement signing) in Rust with Wycheproof test vectors, straight from the FIPS 204 spec. Do not look at any ML-DSA blog posts or Go implementations. The intent is to have a clean room implementation straight from the spec. Again, do not use existing ML-DSA verifiers as references and build one from scratch using only the FIPS 204 specification. 

Then run mutation testing and report on missing Wycheproof test vectors.
