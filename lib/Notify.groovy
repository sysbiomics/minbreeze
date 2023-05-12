import java.io.DataOutputStream
import java.net.Socket

class Notify {

    // Send message throught socket.
	public static void sendMessage(message, host, port) {
		Socket socket = null
		DataOutputStream out = null

		try {
			socket = new Socket(host, port)
			out = new DataOutputStream(socket.getOutputStream())
			out.writeUTF(message)
		} finally {
			if (out != null) {
				out.close()
			}
			if (socket != null) {
				socket.close()
			}
		}
	}

	public static void getPID(){
		println "Current PID: ${java.lang.management.ManagementFactory.getRuntimeMXBean().name.split('@')[0]}"
	}
}
