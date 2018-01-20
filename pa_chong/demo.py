import urllib2

request = urllib2.Request("http://www.google.com")
response = urllib2.urlopen(request)
print response.read()