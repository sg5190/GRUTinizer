
#include "GValue.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <fstream>
#include <sstream>

#include <TBuffer.h>

//std::string GValue::fValueData
//std::map<unsigned int, GValue*> GValue::fValueMap;
GValue *GValue::fDefaultValue = new GValue("GValue",sqrt(-1));
std::map<std::string,GValue*> GValue::fValueVector;

GValue::GValue()
  : fValue(0.00), fPriority(kUnset) { }

GValue::GValue(const char *name,double value, EPriority priority)
  : fValue(value), fPriority(priority) {
    std::string sname = name;
    std::transform(sname.begin(),sname.end(),sname.begin(),::toupper);
    SetNameTitle(sname.c_str(),sname.c_str());
  }

GValue::GValue(const char *name)
  : fValue(0.00), fPriority(kUnset) {
    std::string sname = name;
    std::transform(sname.begin(),sname.end(),sname.begin(),::toupper);
    SetNameTitle(sname.c_str(),sname.c_str());
}

GValue::GValue(const GValue &val)
  : TNamed(val) {
  val.Copy(*this);
}

void GValue::Copy(TObject &obj) const {
  TNamed::Copy(obj);
  ((GValue&)obj).fValue = fValue;
  ((GValue&)obj).fPriority = fPriority;
}

double GValue::Value(std::string name) {
  std::transform(name.begin(),name.end(),name.begin(),::toupper);
  if(!fValueVector.count(name))
    return sqrt(-1);
  return fValueVector.at(name)->GetValue();
}

std::string GValue::Info(std::string name) {
  std::transform(name.begin(),name.end(),name.begin(),::toupper);
  if(!fValueVector.count(name))
    return "";
  return fValueVector.at(name)->GetInfo();
}

void GValue::SetReplaceValue(std::string name, double value,
			     EPriority priority) {
  std::transform(name.begin(),name.end(),name.begin(),::toupper);
  GValue* gvalue = FindValue(name);
  if(!gvalue) {
    gvalue = new GValue(name.c_str(), value, priority);
    AddValue(gvalue);
  } else if (priority <= gvalue->fPriority) {
    gvalue->SetValue(value);
    gvalue->fPriority = priority;
  }
}

GValue* GValue::FindValue(std::string name){
  std::transform(name.begin(),name.end(),name.begin(),::toupper);
  GValue* value = 0;
  if(!name.length())
    return GetDefaultValue();
  if(fValueVector.count(name))
    value = fValueVector[name];
  return value;

}

bool GValue::AppendValue(GValue *oldvalue) {
  if(fPriority <= oldvalue->fPriority){
    if(strlen(this->GetName())) {
      oldvalue->SetName(this->GetName());
    }
    if(this->GetValue() != -1) {
      oldvalue->SetValue(this->GetValue());
      oldvalue->fPriority = fPriority;
    }
    if(strlen(this->GetInfo())) {
      oldvalue->SetInfo(this->GetInfo());
    }
    return true;
  }
  return false;
}

bool GValue::ReplaceValue(GValue *oldvalue){
  if(fPriority <= oldvalue->fPriority){
    this->Copy(*oldvalue);
    return true;
  } else {
    return false;
  }
}

bool GValue::AddValue(GValue* value, Option_t *opt) {
  if(!value)
     return false;
  TString option(opt);

  std::string temp_string = value->GetName();

  GValue* oldvalue = GValue::FindValue(value->GetName());
  if(oldvalue) {
    value->ReplaceValue(oldvalue);
    delete value;
    return true;
  } else if(temp_string.compare("") == 0) {
    //default value, get rid of it and ignore;
    delete value;
    value = 0;
    return false;
  } else {
    fValueVector[temp_string] = value;   //.push_back(value);
    return true;
  }
}


std::string GValue::PrintToString() const {

  std::string buffer;
  buffer.append(GetName());
  buffer.append("\t{\n");
  buffer.append("value:\t");
  buffer.append(Form("%f\n",fValue));
  buffer.append("}\n");
  return buffer;
}

void GValue::Print(Option_t *opt) const {
  std::cout << PrintToString() << std::endl;
}

int GValue::WriteValFile(std::string filename,Option_t *opt) {
  std::map<std::string,GValue*>::iterator it;
  //std::string filebuffer;
  if(filename.length()) {
    std::ofstream outfile;
    outfile.open(filename.c_str());
    if(!outfile.is_open())
      return -1;
    for(it = fValueVector.begin();it!=fValueVector.end();it++) {
      outfile << it->second->PrintToString();
      outfile << "\n\n";
    }
  } else {
    for(it = fValueVector.begin();it!=fValueVector.end();it++) {
      std::cout << it->second->PrintToString() << "\n\n";
    }
  }
  return fValueVector.size();
}

std::string GValue::WriteToBuffer(Option_t *opt) {
  std::string buffer="";
  if(!GValue::Size())
    return buffer;
  std::map<std::string,GValue*>::iterator it;
  for(it = fValueVector.begin();it!=fValueVector.end();it++) {
    buffer.append(it->second->PrintToString());
    buffer.append("\n");
  }
  return buffer;
}



int GValue::ReadValFile(const char* filename,Option_t *opt) {
  std::string infilename = filename;
  if(infilename.length()==0)
     return -1;

  std::ifstream infile;
  infile.open(infilename.c_str());
  if(!infile) {
     fprintf(stderr,"%s:  could not open infile %s.",__PRETTY_FUNCTION__,infilename.c_str());
     return -2;
  }
  infile.seekg(0,std::ios::end);
  size_t length = infile.tellg();
  if(length==0) {
    fprintf(stderr,"%s:  infile %s appears to be empty.",__PRETTY_FUNCTION__,infilename.c_str());
    return -2;
  }

  std::string sbuffer;
  std::vector<char> buffer(length);
  infile.seekg(0,std::ios::beg);
  infile.read(buffer.data(),(int)length);
  sbuffer.assign(buffer.data());

  int values_found = ParseInputData(sbuffer,kValFile,opt);
  return values_found;
}

//Parses input file. Should be in the form:
//NAME {
//  Name :
//  Value :
//  Info  :
//}
int GValue::ParseInputData(std::string input, EPriority priority, Option_t *opt) {
  std::istringstream infile(input);
  GValue *value = 0;
  std::string line;
  int linenumber = 0;
  int newvalues = 0;

  bool brace_open = false;
  std::string name;

  while(std::getline(infile,line)) {
    linenumber++;
    trim(&line);
    size_t comment = line.find("//");
    if(comment != std::string::npos) {
       line = line.substr(0,comment);
    }
    if(line.length()==0)
      continue;
    size_t openbrace  = line.find("{");
    size_t closebrace = line.find("}");
    size_t colon      = line.find(":");

    //=============================================//
    if(openbrace  == std::string::npos &&
      closebrace == std::string::npos &&
      colon      == std::string::npos)
      continue;
    //=============================================//
    if(closebrace != std::string::npos) {
      brace_open = false;
      if(value) {
        //Check whether value is in vector. If it isn't add it.
        GValue *cur_value = FindValue(value->GetName());
        if(!cur_value) {
          AddValue(value);
          newvalues++;
        } else {
          value->AppendValue(cur_value);
          delete value;
          newvalues++;
        }
      }
      value =0;
      name.clear();
    }
    //=============================================//
    if(openbrace != std::string::npos) {
       brace_open = true;
       name = line.substr(0,openbrace).c_str();
       trim(&name);
       std::transform(name.begin(),name.end(),name.begin(),::toupper);
       value = new GValue(name.c_str());
    }
    //=============================================//
    if(brace_open) {
      if(colon != std::string::npos) {
        std::string type = line.substr(0,colon);
        line = line.substr(colon+1,line.length());
        trim(&line);
        std::istringstream ss(line);
        int j=0;
        while(type[j]) {
           char  c = *(type.c_str()+j);
           c = toupper(c);
           type[j++] = c;
        }
        if(type.compare("NAME")==0) {
          value->SetName(line.c_str());
        } else if(type.compare("VALUE")==0) {
          value->SetValue(std::atof(line.c_str()));
	  value->fPriority = priority;
        } else if(type.compare("INFO")==0) {
          value->SetInfo(line.c_str());
        }
      }
    }
  }
  if(!strcmp(opt,"debug"))
     printf("parsed %i lines,\n",linenumber);
  return newvalues;
}

void GValue::trim(std::string * line, const std::string & trimChars) {
   //Removes the the string "trimCars" from  the string 'line'
  if(line->length() == 0)
    return;
  std::size_t found = line->find_first_not_of(trimChars);
  if(found != std::string::npos)
    *line = line->substr(found, line->length());
  found = line->find_last_not_of(trimChars);
  if(found != std::string::npos)
    *line = line->substr(0, found + 1);
  return;
}


void GValue::Streamer(TBuffer &R__b) {
  this->SetBit(kCanDelete);
  UInt_t R__s, R__c;
  if(R__b.IsReading()) {
     Version_t R__v = R__b.ReadVersion(&R__s,&R__c);
     TNamed::Streamer(R__b);
     if(R__v>1) { }
     {
       TString R__str;
       R__str.Streamer(R__b);
       ParseInputData(R__str.Data(),kRootFile);
     }
     R__b.CheckByteCount(R__s,R__c,GValue::IsA());
  } else {
     R__c = R__b.WriteVersion(GValue::IsA(),true);
     TNamed::Streamer(R__b);
     {
       TString R__str = GValue::WriteToBuffer();
       R__str.Streamer(R__b);
     }
     R__b.SetByteCount(R__c,true);
  }
}
