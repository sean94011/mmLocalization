#include <Stepper.h>
#include <Timer.h>

const long stepsPerRevolution = 1600;  // only determines the torque and affect the speed
const int setMyStepper = stepsPerRevolution;
// recommend to set it to 4 or 8 times the Microstep number 
//const int NumPoint = 120;
const int MaxNumPoint = 100;
//const int D = 1248; // 1/2*lambda (1.95mm for 77GHz)
const long D = 400; // 1/2*lambda (2.5mm for 77GHz)
Timer t;  
Stepper myStepper(setMyStepper, 6, 7);
int t_start,t_end;

void setup() {
  myStepper.setSpeed(100); // number of setMyStepper arduino step (setMyStepper/2 actual step) per minute
  Serial.begin(9600); // initialize the serial port:
  
  Serial.println("------------------------------------------");
  Serial.println("-------------- x-plotter (max step 200)-----------------");
  Serial.println("Enter Number of step (1/2*lambda step = 2.5mm)");
  Serial.println("Positive: Move right,  Negative: Move left");
  Serial.println("------------------------------------------");
  Serial.println("Step size D = ");
  Serial.println(D);
}

void loop() {
  while (Serial.available() > 0) {
    int NumStep,pauseInit;
    float pauseStep;
    NumStep = Serial.parseInt();
    pauseStep = Serial.parseFloat();
    pauseInit = Serial.parseInt();
    Serial.print("NumStep = ");
    Serial.println(NumStep);

    t.update();
    
    if (Serial.read() == '\n') {
      if (NumStep > MaxNumPoint){
          Serial.println("Number of step exceeds max number of step.");
          return;
      }
      else{
        delay(pauseInit*1000);
        t_start = millis();
        if(NumStep >= 0){
            for (int k=0; k<NumStep; k++){
              myStepper.step(D); // 1*lambda move left 
              delay(pauseStep*1000);
            }                     
        }
        else{
            for (int k=0; k<-NumStep; k++){
              myStepper.step(-D); // 1*lambda move left 
              delay(pauseStep*1000);
            }                               
        }          
        t_end = millis();
        Serial.println(t_end-t_start);
      }
    }    
  }
  
}

// total time = 4992*60/(20*6400) seconds
