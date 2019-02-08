open Owl_ode.Common
open Owl_ode.Types

type 'a f_t = (float, 'a) M.t -> (float, 'a) M.t -> float -> (float, 'a) M.t

val contact1_damped_s :
  a:(float -> float) ->
  f:'a f_t ->
  dt:float -> 
  (float, 'a) M.t ->
  (float, 'a) M.t ->
  float ->
  (float, 'a) M.t * (float, 'a) M.t * float

val contact2_damped_s :
  a:(float -> float) ->
  f:'a f_t ->
  dt:float -> 
  (float, 'a) M.t ->
  (float, 'a) M.t ->
  float ->
  (float, 'a) M.t * (float, 'a) M.t * float


module S: sig
  type mat = Owl_dense_matrix_s.mat

  module Contact1_damped :
    functor (A : sig val a : float -> float end) -> 
      SolverT with type s = mat * mat
               and type t = mat
               and type output = float array * mat * mat

  module Contact2_damped :
    functor (A : sig val a : float -> float end) ->
      SolverT with type s = mat * mat
               and type t = mat
               and type output = float array * mat * mat

end

module D: sig
  type mat = Owl_dense_matrix_d.mat

  module Contact1_damped :
    functor (A : sig val a : float -> float end) -> 
      SolverT with type s = mat * mat
               and type t = mat
               and type output = float array * mat * mat

  module Contact2_damped :
    functor (A : sig val a : float -> float end) ->
      SolverT with type s = mat * mat
               and type t = mat
               and type output = float array * mat * mat

end