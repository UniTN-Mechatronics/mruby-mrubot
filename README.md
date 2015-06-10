# mruby-mrubot
mruby gem for bsic robotics (more to come)

## `Mrubot::NewtonRaphson`

1-D numerical solver.

```ruby
nr = Mrubot::NewtonRaphson.new
nr.f = lambda {|x| 2*x**2 - 15*x + 5}
nr.fd = lambda {|x| 4*x - 15}
nr.x0 = 100                    # or nr.guess = 100
nr.solve
 #  => 7.1503676271839
nr.solve(-100) {|pos| p pos }
 # [0, -48.180722891566, 21505]
 # [1, -22.32668764496, 5370.4749600813]
 # [2, -9.5100456765923, 1336.8622770656]
 # [3, -3.316013040819, 328.53262269063]
 # [4, -0.60118361215881, 76.732080586048]
 # [5, 0.24574672722334, 14.740597653439]
 # [6, 0.34809249733909, 1.4345819995319]
 # [7, 0.34963202430336, 0.020949313321173]
 # [8, 0.34963237281612, 4.7402865481772e-06]
 #  => 0.34963237281612
```


## `Mrubot::Ballistic`
No-drag ballistics.

```ruby
b = Mrubot::Ballistic.new # => #<Mrubot::Ballistic:0x7fbcd4066b60 @y=0, @x=10, @v=10, @theta=0.78539816339745>
b.max_distance # => 10.193679918451
b.distance(angle) # default to @theta
b.reach_angle(distance) # => angle for reaching distance
b.theta_for_xy
b.theta_for_x(x) # accepts x and assumes @y
b.theta_for_y(y) # accepts y and assumes @x
```

## `ProcessInfo`
- `ProcessInfo.current_mem` gives the current process memory occupation
- `ProcessInfo.peak_mem` gives the current process peak memory occupation
- `ProcessInfo::PID` gives the current process ID
- `ProcessInfo::PPID` gives the current _parent_ process ID