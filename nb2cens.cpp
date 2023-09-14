#include <TMB.hpp>

// from TMB source:
// TMB_BIND_ATOMIC(pbeta, 111, toms708::pbeta(x[0], x[1], x[2], 1, 0) )
TMB_BIND_ATOMIC(pbeta_log,
  111, // need derivatives on all 3 'x' arguments
  atomic::toms708::pbeta(x[0], x[1], x[2], 1, 1)) // lower.tail/log.p = true

template<class Type>
Type pbeta_log(Type q, Type shape1, Type shape2) {
    CppAD::vector<Type> args(4); // last index reserved for derivative order
    args[0] = q;
    args[1] = shape1;
    args[2] = shape2;
    args[3] = 0;
    return pbeta_log(args)[0];
}

template <class Type>
Type pnbinom2_log(Type x, Type size, Type mu) {
  // from R source:
  // pnbinom_mu = pbeta(pr, size, x + 1, lower_tail, log_p);
  // pr = size/(size + mu), 1-pr = mu/(size+mu);
  // check:
  // size <- 1; x <- 3; mu <- 2
  // pnbinom(q = x, size = size, mu = mu, log = TRUE)
  // #> -0.2200619
  // pbeta(size/(size + mu), size, x + 1, log = TRUE)
  // #> -0.2200619
  return pbeta_log(size/(size + mu), size, x + Type(1.));
}

template <class Type>
Type dcensnb2_right(Type x, Type size, Type mu, int give_log = 0) {
  Type ll;
  ll = pnbinom2_log(x - Type(1.), size, mu); // F(lower-1)
  ll = logspace_sub(Type(0.), ll); // 1 - F(lower-1)
  if (give_log)
    return ll;
  else
    return exp(ll);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x);
  DATA_VECTOR(y);
  DATA_IVECTOR(cens);
  PARAMETER(b0);
  PARAMETER(b1);
  PARAMETER(ln_phi);
  Type jnll = 0.;
  for (int i = 0; i < y.size(); i++) {
    if (cens(i)) {
      jnll -= dcensnb2_right(y(i), exp(ln_phi), exp(b0 + b1 * x(i)), true);
    } else {
      Type s1 = b0 + b1 * x(i);
      Type s2 = Type(2.) * s1 - ln_phi;
      jnll -= dnbinom_robust(y(i), s1, s2, true);
    }
  }
  return jnll;
}
