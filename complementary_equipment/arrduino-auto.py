
import select
import socket
import serial
import time

from boltons.socketutils import *

arduino = serial.Serial(port='COM3', baudrate=9600, timeout=90)

def write_read(x):
    arduino.write(bytes(x, 'utf-8'))
    time.sleep(0.05)

s = socket.socket()
s.bind(('', 5311)) #172.26.3.173 #127.0.0.1
print("Socket successfully created !!!")

s.listen()
print("Socket is listening")

c, addr = s.accept()
c = BufferedSocket(c, 90)
print("Got connnection from ", addr)

while(True):
    try: 
        data = c.recv(4096)
        data = data.decode()
        print("Got data from ", addr)
        print(data)
        write_read(data)

    except Timeout:
        print("Timeout in socket, trying to reconnect...\nListening for new client...")
        s.listen()

        c, addr = s.accept()
        c = BufferedSocket(c, 90)
        print("Got connnection from ", addr)
