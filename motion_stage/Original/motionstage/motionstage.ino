#include <AccelStepper.h>

void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
}

String rx_str1 = "";
String rx_str2 = "";
double spd = 0;
double location = 0;
AccelStepper stepper(AccelStepper::FULL4WIRE, 2, 3, 4, 5);

void loop() {

  //input is a csv of spd, displacement
  if(Serial.available() > 0)
  {
    rx_str1 = Serial.readStringUntil(',');
    rx_str2 = Serial.readString();
    spd = toDouble(rx_str1);
    location = toDouble(rx_str2);
  }

  //Currently assume mm/s speed, cm position
  //Multiplier will be different depending on motion stage used 
  //See steps/turn and then the drive length of on turn on linear actuator to define
  //Run to new position will drive motor until it reaches the location given (input is location, not distance to move)
  stepper.setSpeed(round(spd*360));
  stepper.setMaxSpeed(round(spd*360));
  stepper.runToNewPosition(round(location*3600)); 
  delay(20);
}
