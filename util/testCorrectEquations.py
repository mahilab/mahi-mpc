import string



fileheader = ["#pragma once\n",
              "#include <vector>\n",
              "#include <Eigen/Dense>\n",
              "#include <mass_properties.hpp>\n",
              "using namespace mahi::util;\n",
              "\n"
              "inline double sq(double num){\n",
              "\treturn num*num;\n",
              "}\n",
              "\n"]

variables = ["\tconst double q0 = qs[0];\n",
             "\tconst double q1 = qs[1];\n",
             "\tconst double q2 = qs[2];\n",
             "\tconst double q3 = qs[3];\n",
             "\tconst double qd0 = qs[4];\n",
             "\tconst double qd1 = qs[5];\n",
             "\tconst double qd2 = qs[6];\n",
             "\tconst double qd3 = qs[7];\n",
             "\tconst double sintheta0 = sin(q0);\n",
             "\tconst double costheta0 = cos(q0);\n",
             "\tconst double sintheta1 = sin(q1);\n",
             "\tconst double costheta1 = cos(q1);\n",
             "\tconst double sintheta2 = sin(q2);\n",
             "\tconst double costheta2 = cos(q2);\n",
             "\tconst double sintheta3 = sin(q3);\n",
             "\tconst double costheta3 = cos(q3);\n",
             "\n"]

filename = "nathanisdumb"


with open("Equations/"+filename+".txt") as read_file:
    fileContent = read_file.read()
    read_file.close()

mat_name = filename
combinedFile = fileContent.split('[[')
firstPart = combinedFile[0]

firstPart = firstPart.replace("[@","\tconst double var")
firstPart = firstPart.replace(" @","\tconst double var")
firstPart = firstPart.replace("\n@","\tconst double var")
firstPart = firstPart.replace("@","var")

firstPart = firstPart.replace(",",";\n")
firstPart = firstPart.replace("sin(q0)","sintheta0")
firstPart = firstPart.replace("cos(q0)","costheta0")
firstPart = firstPart.replace("sin(q1)","sintheta1")
firstPart = firstPart.replace("cos(q1)","costheta1")
firstPart = firstPart.replace("sin(q2)","sintheta2")
firstPart = firstPart.replace("cos(q2)","costheta2")
firstPart = firstPart.replace("sin(q3)","sintheta3")
firstPart = firstPart.replace("cos(q3)","costheta3")
firstPart = firstPart.replace("="," = ")


secondPart = combinedFile[1]
secondPart = secondPart.replace("@","var")
secondPart = secondPart.replace(",","\n")
secondPart = secondPart.replace("]","")
secondPart = secondPart.replace("[","")
secondPart = secondPart.replace("\n\n","\n")
secondPart = secondPart.replace(" (","(")
secondPart = secondPart.replace(" v","v")

secondPartVec = secondPart.split("\n")
i = 0
secondPart2 = ""
for x in secondPartVec:
    if i < 4:
        secondPartVec[i] = "\t" + mat_name + "(0," + str(i) + ") = " + secondPartVec[i] + ";\n"
    elif i < 8:
        secondPartVec[i] = "\t" + mat_name + "(1," + str(i-4) + ") = " + secondPartVec[i] + ";\n"
    elif i < 12:
        secondPartVec[i] = "\t" + mat_name + "(2," + str(i-8) + ") = " + secondPartVec[i] + ";\n"
    else:
        secondPartVec[i] = "\t" + mat_name + "(3," + str(i-12) + ") = " + secondPartVec[i] + ";\n"
    secondPart2 = secondPart2 + secondPartVec[i]
    i = i+1


# test = content.split('[[')
# print(test[1])
write_file = open("../include/" + mat_name + ".hpp","w")
write_file.writelines(fileheader)

write_file.write("inline Eigen::MatrixXd get_" + mat_name + "(const std::vector<double>& qs){\n")
write_file.write("\tEigen::MatrixXd " + mat_name + " = Eigen::MatrixXd::Zero(4,4); \n\n")

write_file.writelines(variables)

write_file.write(firstPart + "\n")
write_file.write(secondPart2 + "\n")
write_file.write("\treturn " + mat_name + ";\n}")

write_file.close()


print("Hello World")