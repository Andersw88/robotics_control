function x_new=step(inverted_pendulum)

[t,X_new] = rk4(@(t,x) dynamics(inverted_pendulum,t,x),[0 inverted_pendulum.dt_],inverted_pendulum.x_,inverted_pendulum.dt_);

x_new=X_new(end,:)';

