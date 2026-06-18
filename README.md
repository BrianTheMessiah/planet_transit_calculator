# planet_transit_calculator

### What is a Hohmann Transfer Orbit?
- Hohmann Transfer Orbit-Orbital maneuver used to transfer spacecraft between two orbits of different altitudes around a central body 
- Idealized case of the Hohmann Transfer is when the orbits are both circular and *coplanar* 
    - *coplanar* means when they are on the same plane 
- Maneuver is accomplished by placing the craft in elliptical transfer orbit tangential to both the initial and target orbits 
    - Utilizes two burns, one to make the transfer orbit, another to adjust the orbit to match the target 
- Travel times are longer with the maneuver, but less delta-v is required  
- Hohmann-German scientist dude who wrote the Attainability of Celestial Bodies (killer name)
- To travel to a different celestial body, the origin body must be in a certain position relative to the destination body 
    - Space missions must wait for the alignment to occur, opening a launch window 
- ![](images/hohmann_diagram.png)
    - The yellow ellipse diagrams what a Hohmann transfer would look like 
    - The perapsis (the point closest to the central celestial body of the orbit) is where the 1 and 2 ellipses touch
    - The apoapsis (the point farthest from the central body of the orbit) is where the 2 and 3 ellipses touch
    - Remember that the Hohmann transfer can be utilized for burns between two different planets, or it can be utilized to move satellites to different orbits  
- To raise your apoapsis, do a prograde burn at perapsis (moves spacecraft to a higher orbit) 
    - **What does a prograde burn mean?**
        - Prograde burn is when spacecraft fires its engine in the same direction as your current velocity vector 
            - Increases orbital speed 
            - When increasing speed at perapsis, spacecrafts apoapsis gets raised (stretch the far side of the orbit outward) 
            - In Hohmann transfer from lower to higher circular orbit, first burn at periapsis is the prograde burn 
- To lower periapsis, do a retrograde burn at apoapsis (moves you toward a lower orbit) 
    - **What does retrograde burn mean in this case?** 
- Type I and Type II
    - Ideal Hohmann transfer orbit is when a spacecraft can transfer between two orbits at exactly 180 degrees (half an orbit of the Hohmann transfer) 
        - **What does 180 degrees around the primary mean**
        - **What does coplanar mean** 
    - Type I: Orbit traversing less than 180 degrees around the primary  
    - Type II: Orbit traversing more than 180 degrees around the primary 
- **Does the Hohmann transfer assume your already in space?**
    - Hohmann transfer does assume your spacecraft is in space, specifically in a stable circular orbit around the central body 
- **What happens before Hohmann transfer?** 
    - Launch from the surface 
    - Reach Low Earth Orbit (LEO, parking orbit) 
    - Circularize 
    - Perform the Hohmann departure burn at periapsis of that orbit 
    - First burn is one that raises apoapsis to the target orbit 
### Calculation 
    - Total energy of a small body orbiting a much larger body is the sum of its kinetic energy and potential energy 
- This total energy also equals half the potential at the average distance $\alpha$ (semi-major axis)
    - **What is the semi-major axis**
    - $$ E = \frac{mv^2}{2} - \frac{GMm}{r} = - \frac{GMm}{2\alpha} $$
    - E: Total mechanical energy of the orbit (kinetic + potential) 
    - m: Mass of orbiting body (satellite, spacecraft planet) 
    - v: Instantaneous orbital speed at distance r 
    - G: Universal gravitational constant ($$ 6.674 * 10^-11 m^3 kg^{-1} s^{-2} $$)
        - **What is the universal gravitational constant?**
    - M: Mass of the central celestial body (Earth, Sun, Jupiter)
    - r: Instantaneous distance between the two bodies 
    - a: Semi major axis of the orbit 
    - For a circle semi-major axis = radius
    - In an ellipse, semi-major axis is half the long axis  
- **What causes circular orbits?**
    - Gravity naturally wants to pull spacecraft/celestial bodies into an ellipse 
    - Circular orbit is special case where spacecraft speed matches exact value needed for constant-radius freefall
        - Circular orbit are lowest energy, most stable solution of the two-body problem 
    - Circular orbit is simply $$ v = \frac{\sqrt{u}}{r} $$
        - v: Orbital velocity needed for circular orbit at radis r 
        - $$ \mu $$: Standard gravitational parameter of the central celestial body 
            - $$ \mu=GM $$
            - **What is the standard gravitational parameter $$ \mu $$?** 

### Questions 
1. What does $\frac{m^3}{s^2}$ mean
2. 
