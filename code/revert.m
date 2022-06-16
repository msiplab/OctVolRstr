SAIVDR_VERSION = "4.2.2.1";
SAIVDR_ROOT = "SaivDr-"+SAIVDR_VERSION;
if exist("./"+SAIVDR_ROOT,"dir")
    rmdir("./"+SAIVDR_ROOT,"s")
end
if exist("../data/materials","dir")
    rmdir("../data/materials","s")
end