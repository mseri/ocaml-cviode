let damped_noforcing xs _ =
  Owl.Mat.mul_scalar xs (-1.0)

let damped_forcing beta omega xs t = 
  Owl.Mat.(add_scalar
             (damped_noforcing xs 0.0)
             (beta *. Owl.Maths.sin(omega *. t))
          )

let a _ = 1.0 
let dt = 0.1

let plot_sol fname t sol1 sol2 =
  let open Owl in
  let h = Plot.create fname in
  let open Plot in
  set_foreground_color h 0 0 0;
  set_background_color h 255 255 255;
  plot ~h ~spec:[ RGB (0,0,255); LineStyle 1 ] t (Mat.col sol1 0);
  plot ~h ~spec:[ RGB (0,255,0); LineStyle 1 ] t (Mat.col sol2 0);
  (* XXX: I could not figure out how to make the legend black instead of red *)
  legend_on h ~position:NorthEast [|"Contact 1st"; "Contact 2nd";|];
  output h

let () =
  let y0 = Owl.Mat.of_array [|-0.25; 0.75|] 1 2 in
  let tspan = (0.0, 15.0) in
  let t = Owl.Arr.linspace 0.0 15.0 (int_of_float @@ Float.floor (15.0/.dt)) in
  let sol1 = Cviode.contact1_damped ~a ~f:damped_noforcing y0 tspan dt in
  let sol2 = Cviode.contact2_damped ~a ~f:damped_noforcing y0 tspan dt in
  plot_sol "damped.png" t sol1 sol2;

  let y0 = Owl.Mat.of_array [|-0.284; -0.027|] 1 2 in
  let tspan = (0.0, 50.0) in
  let t = Owl.Arr.linspace 0.0 50.0 (int_of_float @@ Float.floor (50.0/.dt)) in
  let damped_forcing = damped_forcing 0.3 (Float.pi/.3.0) in
  let sol1 = Cviode.contact2_damped ~a ~f:damped_forcing y0 tspan dt in
  let sol2 = Cviode.contact2_damped ~a ~f:damped_forcing y0 tspan dt in
  plot_sol "forced.png" t sol1 sol2