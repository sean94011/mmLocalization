https://www.arduino.cc/reference/en/libraries/stepper/

Microstep number eg. 3200 on the Microstep Driver
times 4 equals to the actual steps need for a 360 rotation
but in arduino, the step funciton input has to be double of the actual steps

stepsPerRevolution determines the torque, fewer stepsPerRevolution gives more torque, recommend to be 4 times the Microstep number 
 
The setSpeed input v times 2 equals round per minute, but maybe slightly off with different Microstep setting. 
Recommend find out the actual value with Arduino timer library.
  