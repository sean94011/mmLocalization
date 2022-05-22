#include <Stepper.h>

const long stepsPerRevolution=800*8;  // only determines the torque and affect the speed
const int setMyStepper = stepsPerRevolution;
// recommend to set it to 4 or 8 times the Microstep number 
//const int NumPoint = 120;
const int MaxNumPoint = 200;
//const int D = 1248; // 1/2*lambda (1.95mm for 77GHz)
const int D = 1600; // 1/2*lambda (2.5mm for 77GHz)
  
Stepper myStepper(setMyStepper, 6, 7);

void setup() {
  myStepper.setSpeed(75); // number of setMyStepper arduino step (setMyStepper/2 actual step) per minute
  Serial.begin(9600); // initialize the serial port:
  
  Serial.println("------------------------------------------");
  Serial.println("-------------- x-plotter (max step 200)-----------------");
  Serial.println("Enter Number of step (1/2*lambda step = 2.5mm)");
  Serial.println("Positive: Move right,  Negative: Move left");
  Serial.println("------------------------------------------");

}

void loop() {
  while (Serial.available() > 0) {
    int NumStep;
    NumStep = Serial.parseInt();
    Serial.print("NumStep = ");
    Serial.println(NumStep);
    
    if (Serial.read() == '\n') {
      if (NumStep > MaxNumPoint){
          Serial.println("Number of step exceeds max number of step.");
          return;
      }
      else{
        if(NumStep >= 0){
            for (int k=0; k<NumStep; k++){
              myStepper.step(D); // 1*lambda move left 
            }                     
        }
        else{
            for (int k=0; k<-NumStep; k++){
              myStepper.step(-D); // 1*lambda move left 
            }                               
        }          
      }
    }    
  }
  
}

// total time = 4992*60/(20*6400) seconds
