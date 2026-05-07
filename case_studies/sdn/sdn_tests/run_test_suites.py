import subprocess
import shlex
import json
import requests 
import time

device1 = "of:0000000000000001"
device2 = "of:0000000000000002"
app = 666

pattern = ''' curl -u onos:rocks -X POST --header 'Content-Type: application/json' --header 'Accept: application/json' -d '''
url = ''' 'http://localhost:8181/onos/v1/flows?appId=666' '''


payload_for_test = """{
        "flows": [
        {
          "priority": 40000,
          "timeout": 0,
          "isPermanent": false,
          "deviceId": "of:0000000000000001",
          "treatment": {
            "instructions": [
              {
                "type": "OUTPUT",
                "port": "2"
              }
            ]
          },
          "selector": {
            "criteria": [
              {
                "type": "IN_PORT",
                "port": "1"
              },
              {
                "type": "ETH_SRC",
                "mac": "9a:d8:73:d8:90:6a"
              },
              {
                "type": "ETH_DST",
                "mac": "9a:d8:73:d8:90:6b"
              }
            ]
          }
        }
        ]
    }"""

class ONOS_interface:
    def __init__(self):
        self.flow_rules = list()
        return

    def run_curl_cmd(self, cmd):
        args = shlex.split(cmd)
        process = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        resp = json.loads(stdout.decode('utf-8'))
        return resp

    def delete_rules_from_flow_rules(self):
        delete_str = 'http://localhost:8181/onos/v1/flows/of%3A0000000000000001/'
        for flowId in self.flow_rules:
            delete_str_flowId = delete_str + flowId
            resp = requests.delete(delete_str_flowId, auth=('onos', 'rocks'))
            #print("resp =", resp)
        return

    def get_rules_from_flow_rules(self):
        get_str = 'http://localhost:8181/onos/v1/flows/of%3A0000000000000001/'
        for flowId in self.flow_rules:
            get_str_flowId = get_str + flowId
            resp = requests.get(get_str_flowId, auth=('onos', 'rocks'))
            #print("resp =", resp)
        return

    def add_rules(self, rules_number):
        dummy_payload = """{
            "flows": [
            {
              "priority": 40000,
              "timeout": 0,
              "isPermanent": true,
              "deviceId": "of:0000000000000001",
              "treatment": {
                "instructions": [
                  {
                    "type": "OUTPUT",
                    "port": "2"
                  }
                ]
              },
              "selector": {
                "criteria": [
                  {
                    "type": "IN_PORT",
                    "port": "1"
                  },
                  {
                    "type": "ETH_SRC",
                    "mac": "9a:d8:73:d8:90:6a"
                  },
                  {
                    "type": "ETH_DST",
                    "mac": "9a:d8:73:d8:90:6b"
                  }
                ]
              }
            }
            ]
        }"""
        #time.sleep(2)
        data = json.loads(dummy_payload)
        for i in range(0, rules_number):
            data["flows"][0]['priority'] = 40000 + (i+1)
            json_data_payload = json.dumps(data)
            cmd_dummy_payload = pattern + ''' ' ''' + json_data_payload + ''' ' ''' + str(url)
            resp = self.run_curl_cmd(cmd_dummy_payload)
            self.flow_rules.append(resp['flows'][0]['flowId'])
        return

    def run_test_pf1_pf2(self, t1, timeout1, t2, timeout2):
        dummy_payload = """{
            "flows": [
            {
              "priority": 40000,
              "timeout": 0,
              "isPermanent": false,
              "deviceId": "of:0000000000000001",
              "treatment": {
                "instructions": [
                  {
                    "type": "OUTPUT",
                    "port": "2"
                  }
                ]
              },
              "selector": {
                "criteria": [
                  {
                    "type": "IN_PORT",
                    "port": "1"
                  },
                  {
                    "type": "ETH_SRC",
                    "mac": "9a:d8:73:d8:90:6a"
                  },
                  {
                    "type": "ETH_DST",
                    "mac": "9a:d8:73:d8:90:6b"
                  }
                ]
              }
            }
            ]
        }"""
        data1 = json.loads(dummy_payload)
        data1["flows"][0]['priority'] = 40000 + 1
        data1["flows"][0]['timeout'] = int(timeout1)
        json_data_payload1 = json.dumps(data1)
        cmd_dummy_payload1 = pattern + ''' ' ''' + json_data_payload1 + ''' ' ''' + str(url)
        data2 = json.loads(dummy_payload)
        data2["flows"][0]['priority'] = 40000 + 2
        data2["flows"][0]['timeout'] = int(timeout2)
        json_data_payload2 = json.dumps(data2)
        cmd_dummy_payload2 = pattern + ''' ' ''' + json_data_payload2 + ''' ' ''' + str(url)
        # start_test
        time.sleep(t1)
        resp1 = self.run_curl_cmd(cmd_dummy_payload1)
        time.sleep(t2)
        resp2 = self.run_curl_cmd(cmd_dummy_payload2)
        print("resp1=", resp1)
        print("resp2=", resp2)
        return

    def run_pf1(self, timeout1=5):
        data1 = json.loads(payload_for_test)
        data1["flows"][0]['priority'] = 40001
        data1["flows"][0]['timeout'] = int(timeout1)
        json_data_payload1 = json.dumps(data1)
        cmd_dummy_payload1 = pattern + ''' ' ''' + json_data_payload1 + ''' ' ''' + str(url)
        resp1 = self.run_curl_cmd(cmd_dummy_payload1)
        return resp1

    def run_pf2(self, timeout2=2):
        data2 = json.loads(payload_for_test)
        data2["flows"][0]['priority'] = 40002
        data2["flows"][0]['timeout'] = int(timeout2)
        json_data_payload2 = json.dumps(data2)
        cmd_dummy_payload2 = pattern + ''' ' ''' + json_data_payload2 + ''' ' ''' + str(url)
        resp2 = self.run_curl_cmd(cmd_dummy_payload2)
        return resp2

    def run_pf(self, priority, timeout):
        data = json.loads(payload_for_test)
        data["flows"][0]['priority'] = int(priority)
        data["flows"][0]['timeout'] = int(timeout)
        json_data_payload = json.dumps(data)
        cmd_dummy_payload = pattern + ''' ' ''' + json_data_payload + ''' ' ''' + str(url)
        resp = self.run_curl_cmd(cmd_dummy_payload)
        return resp

class TTT:
    def __init__(self, onos):
        self.onos = onos
        return

    def run_seq_1(self, t1, t2, t3, t4):
        self.onos.run_pf1()
        time.sleep(t2-t1)
        self.onos.run_pf1()
        time.sleep(t3-t2)
        self.onos.run_pf2()
        time.sleep(t4-t3)
        self.onos.run_pf1()
        return

    def run_seq_2(self, t1, t2, t3, t4):
        self.onos.run_pf2()
        time.sleep(t2-t1)
        self.onos.run_pf1()
        time.sleep(t3-t2)
        self.onos.run_pf2()
        time.sleep(t4-t3)
        self.onos.run_pf2()
        return

    def run_ttt_mean(self):
        self.run_seq_1(8.0, 18.0, 20.0, 28.5)
        time.sleep(10.0)
        self.run_seq_2(8.0, 9.5, 19.5, 28.0)
        time.sleep(10.0)
        return

    def run_ttt_random(self):
        self.run_seq_1(3.83, 16.52, 17.94, 31.29)
        time.sleep(10.0)
        self.run_seq_2(5.90, 7.78, 14.20, 17.05)
        time.sleep(10.0)
        return

class UTS:
    def __init__(self, onos):
        self.onos = onos
        return

    def run_seq_1(self, t1, t2, t3, t4):
        self.onos.run_pf1()
        time.sleep(t2-t1)
        self.onos.run_pf1()
        time.sleep(t3-t2)
        self.onos.run_pf2()
        time.sleep(t4-t3)
        self.onos.run_pf1()
        return

    def run_seq_2(self, t1, t2, t3, t4):
        self.onos.run_pf2()
        time.sleep(t2-t1)
        return

class TEST_SUITE:
    def __init__(self, onos, file_name):
        self.onos = onos
        self.file = open(file_name, 'r')
        return

    def run_test_sute():
        return
    

#print("(", i, end=",")
#print(t, end=")")
#print()

def run_test_sute(onos, file_name, test_type):
    file = open(file_name, 'r')
    line_number = 0
    for line in file:
        line_list = line.strip().split()
        j = 0
        curr_t = 0
        for i_t in line_list:
            i_s = i_t.strip().split(',')[0]
            t_s = i_t.strip().split(',')[1]
            i = int(i_s[1:])
            t = float(t_s[:len(t_s)-1])
            if j:
                time.sleep(t - curr_t)
            if test_type == "ttt_mean":
                if i == 0:
                    onos.run_pf(1000, 10)
                else:
                    onos.run_pf(1001, 7)
            elif test_type == "ttt_random":
                if i == 0:
                    onos.run_pf(2000, 10)
                else:
                    onos.run_pf(2001, 7)
            elif test_type == "uts_mean":
                if i == 0:
                    onos.run_pf(3000, 10)
                else:
                    onos.run_pf(3001, 7)
            elif test_type == "uts_random":
                if i == 0:
                    onos.run_pf(4000, 10)
                else:
                    onos.run_pf(4001, 7)
            elif test_type == "uppaal_edges":
                if i == 0:
                    onos.run_pf(5000, 10)
                else:
                    onos.run_pf(5001, 7)
            elif test_type == "uppaal_locations":
                if i == 0:
                    onos.run_pf(6000, 10)
                else:
                    onos.run_pf(6001, 7)
            curr_t = float(t)
            j += 1
        time.sleep(10.0)
        print("line_number =", line_number)
        line_number += 1
    return

onos = ONOS_interface()

print("ttt_mean:")
run_test_sute(onos, "shorts_tests/ttt_mean.txt", "ttt_mean")
time.sleep(5.0)
print("ttt_random:")
run_test_sute(onos, "shorts_tests/ttt_random.txt", "ttt_random")
time.sleep(5.0)
print("uts_mean:")
run_test_sute(onos, "shorts_tests/uts_mean.txt", "uts_mean")
time.sleep(5.0)
print("uts_random:")
run_test_sute(onos, "shorts_tests/uts_random.txt", "uts_random")
time.sleep(5.0)
print("uppaal_edges:")
run_test_sute(onos, "shorts_tests/uppaal_edges.txt", "uppaal_edges")
time.sleep(5.0)
print("uppaal_locations:")
run_test_sute(onos, "shorts_tests/uppaal_locations.txt", "uppaal_locations")
time.sleep(5.0)